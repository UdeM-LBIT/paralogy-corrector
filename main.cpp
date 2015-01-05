#define NO_ENSEMBL

/*#include <QtCore/QCoreApplication>*/
#include "trees/paralogycorrector.h"

#include "trees/newicklex.h"
#include "div/util.h"

#include <string>
#include <vector>
/*#include <QDebug>
#include <QProcess>*/

#include <iostream>
#include <fstream>

using namespace std;




/**
 What's going on here :
 Every tree we want to correct is in a .parc file (here, the hard-coded test.parc),
 whcih is a pseudo-xml file.
 This file includes, for every tree, the list of gene names,
 their corresponding species, the gene tree, the species tree,
 and the orthology constraints.  Most of the main() is dedicated to parsing this file,
 and the magic happens on these two lines :
    ParalogyCorrector pc;
    ctree = pc.CorrectGeneTree(genetree, speciestree, geneSpeciesMapping, orthologs);
 The .parc file is not necessary, but here it is used to fill up
 the genetree, speciestree, geneSpeciesMapping, orthologs parameters.
 */
int main(int argc, char *argv[])
{



    //leave this set empty to parse all trees
    set<string> restrictToTrees;


    map<string, string> args;
    string format = "new";

    string prevArg = "";
    for (int i = 0; i < argc; i++)
    {
        if (prevArg != "" && prevArg[0] == '-')
        {
            args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
        }

        prevArg = string(argv[i]);

    }


    if (args.find("i") == args.end())
    {
        cout<<"Error : please specify an input .parc file with -i"<<endl;
        return -1;
    }

    string infile = args["i"];

    string outfile = "";
    if (args.find("o") != args.end())
    {
        outfile = args["o"];
    }

    if (args.find("f") != args.end())
    {
        format = args["f"];
    }

    if (!Util::FileExists(infile))
    {
        cout<<"Error : specified infile "<<infile<<" does not exist"<<endl;
        return -1;
    }

    std::ifstream infilestream(infile);

    //Parse the parc file, extract everythng
    //allModeStrings is the list of instances to correct (demarcated by <INSTANCE>...</INSTANCE>)
    //The map<string, string> for a given instance in allModeStrings,
    //contains specific values of the instance, the key being the xml tag and the value, well, the value
    //eg. string newickk = my_map["<GENETREE>"]
    vector<map<string, string> > allModeStrings;
    map<string, string> curModeStrings;

    bool inInstance = false;
    std::string line;
    string mode = "";
    string modestr = "";
    while (std::getline(infilestream, line))
    {
        if (mode == "")
        {
            if (line == "<INSTANCE>")
                inInstance = true;
            else if (line == "</INSTANCE>" && inInstance)
            {

                if (restrictToTrees.size() == 0 || restrictToTrees.find(curModeStrings["<TREEID>"]) != restrictToTrees.end())
                {
                    allModeStrings.push_back(curModeStrings);
                }

                curModeStrings = map<string, string>();
                mode = "";
                modestr = "";
                inInstance = false;
            }
            else if (line == "<CONSTRAINTS>" || line == "<GENETREE>" || line == "<SPECIESTREE>" || line == "<GENESPECIESMAPPING>" || line == "<TREEID>")
                mode = line;
        }
        else
        {

            string ender = Util::ReplaceAll(mode, "<", "</");
            if (line == ender)
            {
                curModeStrings[mode] = modestr;
                mode = "";
                modestr = "";
            }
            else
            {
                if (line != "")
                {
                    if (modestr != "")
                        modestr += "\n";
                    modestr += line;
                }
            }

        }
    }

    infilestream.close();


    string strout = "";

    //For each <INSTANCE> in the .parc file, build up the parameters and correct the tree
    for (int i = 0; i < allModeStrings.size(); i++)
    {
        map<string, string> modeStrings = allModeStrings[i];

        string treeID = modeStrings["<TREEID>"];
        Node* genetree = NewickLex::ParseNewickString(modeStrings["<GENETREE>"], false);
        Node* speciestree = NewickLex::ParseNewickString(modeStrings["<SPECIESTREE>"], true);
        Node* ctree = NULL;


        //Here we parse the <GENESPECIESMAPPING> to build, you guessed it, geneSpeciesMapping
        map<Node*, Node*> geneSpeciesMapping;
        vector<string> geneSpeciesLines = Util::Split(modeStrings["<GENESPECIESMAPPING>"], "\n");

        //TODO : this lazy loop really sucks - we re-search the tree everytime we need a specific node
        //       Consider using treeInfo, which keeps a map of label -> node
        for (int j = 0; j < geneSpeciesLines.size(); j++)
        {
            vector<string> gs = Util::Split(geneSpeciesLines[j], ":");
            Node* g = genetree->GetNodeWithLabel(gs[0], true);
            Node* s = speciestree->GetNodeWithLabel(gs[1], true);

            geneSpeciesMapping[g] = s;
        }


        //And here we parse the "<CONSTRAINTS>", the pairs of gene names required to be orthologs
        vector<pair<string, string> > orthologs;
        vector<string> orthoLines = Util::Split(modeStrings["<CONSTRAINTS>"], "\n");

        for (int j = 0; j < orthoLines.size(); j++)
        {
            vector<string> xz = Util::Split(orthoLines[j], ":");
            pair<string, string> p = make_pair(xz[0], xz[1]);
            orthologs.push_back(p);
        }


        ParalogyCorrector pc;
        ctree = pc.CorrectGeneTree(genetree, speciestree, geneSpeciesMapping, orthologs);

        //ctree might be NULL if no correction was applied
        if (ctree)
        {
            if (format == "old")
            {
                strout += "TREEID : " + treeID + "\n";
                strout += "BEFORE : " + NewickLex::ToNewickString(genetree) + "\n";
                strout += "AFTER : " + NewickLex::ToNewickString(ctree) + "\n";
            }
            else
            {
                strout += "<INSTANCE>\n<TREEID>\n" + treeID + "\n</TREEID>\n<BEFORE>\n" + NewickLex::ToNewickString(genetree) +
                        "\n</BEFORE>\n<AFTER>\n" + NewickLex::ToNewickString(ctree) + "\n</AFTER>\n</INSTANCE>\n";

            }

            delete ctree;
        }
        else
        {
            if (format == "old")
            {
                strout += "TREEID : " + treeID + " WAS NOT CORRECTED";
            }
            else
            {
                strout += "<INSTANCE>\n<TREEID>\n" + treeID + "\n</TREEID><BEFORE>\n" + NewickLex::ToNewickString(genetree) +
                        "\n</BEFORE>\n<AFTER>\n" + NewickLex::ToNewickString(genetree) + "\n</AFTER>\n<ERROR>\nWas not corrected\n</ERROR>\n</INSTANCE>\n";

            }
        }


        geneSpeciesMapping.clear();
        geneSpeciesLines.clear();
        orthologs.clear();
        orthoLines.clear();

        if (outfile != "")
        {
            ofstream outfilestream;
            outfilestream.open (outfile);
            outfilestream << strout;
            outfilestream.close();
        }
        else
        {
            cout<<strout;
        }

        delete genetree;
        delete speciestree;



    }



    return 0;
}





//What follows is just the implementation of the pipeline used to correct a bunch of trees and evaluate them
//BEWARE THE HARD-CODED PATHS
/*string GetFastaFetcherCommand(Node* tree, string &outfile)
{
    //WE ASSUME GENE LABEL FORMAT IS species_geneid

    //the getfasta command invokes a perl script that fetches all nucleotide sequences of the tree genes
    //this command will be dropped if user chose to use peptide sequences instead
    string commandFasta = "getfasta_b.perl -s \"";

    TreeIterator* it = tree->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), "_");

            commandFasta += sg[0] + ":" + sg[1] + ";";
            //commandFasta += EnsemblUtil::Instance()->GetSpeciesEnsemblName(info->taxa, GLOB_COMPARA_DBNAME) + ":" + info->geneStableID + ";";

        }
    }
    tree->CloseIterator(it);


    commandFasta += "\"";

    if (outfile != "")
        commandFasta += " > " + outfile;

    commandFasta+= "\n";

    //qDebug()<<QString::fromStdString(commandFasta);

    return commandFasta;

}


void FetchFasta(Node* genetree, string &outfile, string &treeID)
{

    //fasta file might exist
    ifstream ifile(outfile);
    if (ifile) {
        qDebug()<<QString::fromStdString(outfile)<<" exists";
        ifile.close();
        return;
    }


    //fasta file might exist from another project
    string prevFasta = "/u/lafonman/tmp/temp" + treeID + ".fasta";
    ifstream tfile(prevFasta, ios::binary);
    if (tfile) {
        qDebug()<<QString::fromStdString(prevFasta)<<" exists (will copy)";


        std::ofstream ofs(outfile, std::ios::binary);

        ofs << tfile.rdbuf();

        ofs.close();
        tfile.close();
        return;
    }



    string tstr = "";
    string command = GetFastaFetcherCommand(genetree, tstr);

    QProcess *fastaProcess = new QProcess();

    fastaProcess->setStandardOutputFile(QString::fromStdString(outfile));

    qDebug()<<"FASTA-ING" <<QString::fromStdString(treeID + " -- " + command);

    fastaProcess->start(QString::fromStdString(command));
    fastaProcess->waitForFinished(-1);


    qDebug()<<"DONE FASTA-ING WITH CODE "<<fastaProcess->error();

    delete fastaProcess;
}


void AlignFile(string &fastaFile, string &outfile, string &repairedOutfile)
{
    QProcess *clustalProcess = new QProcess();


    string command = "clustalw2 -INFILE=" + fastaFile + " -OUTFILE=" + outfile + " -OUTPUT=NEXUS";
    qDebug()<<"ALIGNING "; //<<QString::fromStdString(command);

    clustalProcess->start(QString::fromStdString(command));
    clustalProcess->waitForFinished(-1);


    qDebug()<<"DONE ALIGNING WITH CODE "<<clustalProcess->error();

    delete clustalProcess;






    QProcess *repairProcess = new QProcess();

    string in = outfile;


    command = "nexrepair.py -i \"" + in + "\" -o \"" + repairedOutfile + "\"";
    qDebug()<<"REPAIRING";

    repairProcess->start(QString::fromStdString(command));
    repairProcess->waitForFinished(-1);


    qDebug()<<"DONE REPAIRING WITH CODE "<<repairProcess->error();

    delete repairProcess;

}





void ExecutePhyML(string alignedFile, string treesFile)
{
    ifstream alfile(alignedFile);
    if (!alfile) {
        qDebug()<<QString::fromStdString(alignedFile)<<" DOESN'T exists FOR phyml !";
        return;
    }
    alfile.close();

    string command = "phyml -i " + alignedFile +
                     " -u " + treesFile + " -n 1 -o l --print_site_lnl --no_memory_check";

    QProcess *phymlProcess = new QProcess();


    qDebug()<<"PHYML-ing "<<QString::fromStdString(command);


    phymlProcess->setProcessChannelMode(QProcess::MergedChannels);

    phymlProcess->start(QString::fromStdString(command));


    phymlProcess->waitForReadyRead();
    while (phymlProcess->canReadLine())
    {
        QByteArray bline = phymlProcess->readLine();
        QString line(bline);

        if (line.toUpper().contains("TYPE ANY KEY"))
        {
            phymlProcess->putChar('13');
        }
    }
    phymlProcess->waitForFinished(-1);


    qDebug()<<"DONE PHYML-ing WITH CODE "<<phymlProcess->error();

    delete phymlProcess;
}



void ExecuteConsel(string &phyMLFile, string &resultsFile)
{
    string command = "/u/lafonman/consel/consel/bin/makermt --phyml " + phyMLFile;

    QProcess *mkmtProcess = new QProcess();


    qDebug()<<"MAKERMT";

    mkmtProcess->start(QString::fromStdString(command));
    mkmtProcess->waitForFinished(-1);


    qDebug()<<"DONE MAKERMT WITH CODE "<<mkmtProcess->error();

    delete mkmtProcess;




    string noext = Util::ReplaceAll(phyMLFile, ".txt", "");
    command = "/u/lafonman/consel/consel/bin/consel " + noext;

    QProcess *conselProcess = new QProcess();


    qDebug()<<"CONSEL";

    conselProcess->start(QString::fromStdString(command));
    conselProcess->waitForFinished(-1);


    qDebug()<<"DONE MAKERMT WITH CODE "<<conselProcess->error();

    delete conselProcess;



    command = "/u/lafonman/consel/consel/bin/catpv " + noext + ".pv -s 1";

    QProcess *catpvProcess = new QProcess();


    qDebug()<<"CATPV";

    catpvProcess->setStandardOutputFile(QString::fromStdString(resultsFile));
    catpvProcess->start(QString::fromStdString(command));
    catpvProcess->waitForFinished(-1);



    qDebug()<<"DONE CATPV WITH CODE "<<catpvProcess->error();

    delete catpvProcess;

}







void OutputTrees(Node* originalTree, Node* correctedTree, string &outfile, string &treeID)
{
    //WE ASSUME GENE LABEL FORMAT IS species_geneid

    //make a copy since we change labels
    Node* treeBefore = new Node(false);
    treeBefore->CopyFrom(originalTree);

    Node* treeAfter = new Node(false);
    treeAfter->CopyFrom(correctedTree);

    TreeIterator* it = treeBefore->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), "_");
            n->SetLabel(sg[1]);
        }
        else
        {
            n->SetLabel("");
        }
    }
    treeBefore->CloseIterator(it);

    it = treeAfter->GetPostOrderIterator();
    while (Node* n = it->next())
    {
        if (n->IsLeaf())
        {
            vector<string> sg = Util::Split(n->GetLabel(), "_");
            n->SetLabel(sg[1]);
        }
        else
        {
            n->SetLabel("");
        }
    }
    treeBefore->CloseIterator(it);


    ofstream myfile;
    myfile.open (outfile);
    myfile << NewickLex::ToNewickString(treeBefore) << "\n" << NewickLex::ToNewickString(treeAfter);
    myfile.close();



    delete treeBefore;
    delete treeAfter;

}


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);


    //leave this set empty to parse all trees
    set<string> restrictToTrees;

    //when outfasta != "", it just outputs the fasta fetch command for the given tree only, and does nothing else
    string outfasta = ""; //ENSGT00390000008140";
    string outfastaFile = "/u/lafonman/Projects/ParalogyCorrector/work/sequences/" + outfasta + ".fasta";


    string treesDir = "/u/lafonman/Projects/ParalogyCorrector/work/trees/";
    string fastaDir = "/u/lafonman/Projects/ParalogyCorrector/work/unaligned_seqs/";
    string alignedDir = "/u/lafonman/Projects/ParalogyCorrector/work/aligned_seqs/";
    string conselOutDir = "/u/lafonman/Projects/ParalogyCorrector/work/results/";


    std::ifstream infile("/u/lafonman/Projects/ParalogyCorrector/work/v2-3.parc");



    //Parse the parc file, extract everythng
    vector<map<string, string> > allModeStrings;
    map<string, string> curModeStrings;

    bool inInstance = false;
    std::string line;
    string mode = "";
    string modestr = "";
    while (std::getline(infile, line))
    {
        if (mode == "")
        {
            if (line == "<INSTANCE>")
                inInstance = true;
            else if (line == "</INSTANCE>" && inInstance)
            {
                if (outfasta == "" || outfasta == curModeStrings["<TREEID>"])
                {
                    if (restrictToTrees.size() == 0 || restrictToTrees.find(curModeStrings["<TREEID>"]) != restrictToTrees.end())
                    {
                        allModeStrings.push_back(curModeStrings);
                    }
                }
                curModeStrings = map<string, string>();
                mode = "";
                modestr = "";
                inInstance = false;
            }
            else if (line == "<CONSTRAINTS>" || line == "<GENETREE>" || line == "<SPECIESTREE>" || line == "<GENESPECIESMAPPING>" || line == "<TREEID>")
                mode = line;
        }
        else
        {

            string ender = Util::ReplaceAll(mode, "<", "</");
            if (line == ender)
            {
                curModeStrings[mode] = modestr;
                mode = "";
                modestr = "";
            }
            else
            {
                if (line != "")
                {
                    if (modestr != "")
                        modestr += "\n";
                    modestr += line;
                }
            }

        }
    }

    infile.close();



    string outstr;
    QString qoutstr;
    int nbNull = 0;

    for (int i = 0; i < allModeStrings.size(); i++)
    {
        map<string, string> modeStrings = allModeStrings[i];

        string treeID = modeStrings["<TREEID>"];
        Node* genetree = NewickLex::ParseNewickString(modeStrings["<GENETREE>"], false);
        Node* speciestree = NewickLex::ParseNewickString(modeStrings["<SPECIESTREE>"], true);
        Node* ctree = NULL;

        //special mode - just get fasta sequence
        if (outfasta != "")
        {
            GetFastaFetcherCommand(genetree, outfastaFile);
        }
        else
        {
            map<Node*, Node*> geneSpeciesMapping;
            vector<string> geneSpeciesLines = Util::Split(modeStrings["<GENESPECIESMAPPING>"], "\n");

            //TODO : this lazy loop really sucks - we re-search the tree everytime we need a specific node
            //       Consider using treeInfo, which keeps a map of label -> node
            for (int j = 0; j < geneSpeciesLines.size(); j++)
            {
                vector<string> gs = Util::Split(geneSpeciesLines[j], ":");
                Node* g = genetree->GetNodeWithLabel(gs[0], true);
                Node* s = speciestree->GetNodeWithLabel(gs[1], true);

                geneSpeciesMapping[g] = s;
            }


            vector<pair<string, string> > orthologs;
            vector<string> orthoLines = Util::Split(modeStrings["<CONSTRAINTS>"], "\n");

            for (int j = 0; j < orthoLines.size(); j++)
            {
                vector<string> xz = Util::Split(orthoLines[j], ":");
                pair<string, string> p = make_pair(xz[0], xz[1]);
                orthologs.push_back(p);
            }

            string conselOut = conselOutDir + treeID + ".consel";
            ifstream resultsfile(conselOut);
            if (resultsfile) {
                qDebug()<<"RESULTS " + QString::fromStdString(conselOut) + " exists";
                resultsfile.close();
            }
            else
            {

                ParalogyCorrector pc;
                ctree = pc.CorrectGeneTree(genetree, speciestree, geneSpeciesMapping, orthologs);


                if (ctree)
                {
                    //qDebug()<<QString::number(i) + "/" + QString::number(allModeStrings.size());

                    //string newick = NewickLex::ToNewickString(ctree);

                    //for some extremely weird reasons, that line below makes the program enter an infinite loop
                    //outstr += treeID + "\n" + newick + "\n\n";

                    //qoutstr += QString::fromStdString(treeID + "\n" + newick + "\n\n");

                    if (treesDir != "")
                    {
                        qDebug()<<"Outputting "<<QString::fromStdString(treesDir + treeID)<<".tree";



                        string treesFile = treesDir + treeID + ".trees";
                        OutputTrees(genetree, ctree, treesFile, treeID);

                        if (fastaDir != "")
                        {
                            string fastaFile = fastaDir + treeID + ".fasta";
                            FetchFasta(genetree, fastaFile, treeID);

                            if (alignedDir != "")
                            {
                                string aligned = alignedDir + treeID + ".nxs";
                                string repaired = alignedDir + treeID + "-rep.nxs";

                                AlignFile(fastaFile, aligned, repaired);


                                ExecutePhyML(repaired, treesFile);



                                if (conselOutDir != "")
                                {
                                    string conselIn = repaired + "_phyml_lk.txt";

                                    ExecuteConsel(conselIn, conselOut);
                                }

                            }

                        }
                    }

                    delete ctree;
                }
                else
                {
                    nbNull++;
                }

            }


            geneSpeciesMapping.clear();
            geneSpeciesLines.clear();
            orthologs.clear();
            orthoLines.clear();
        }

        delete genetree;
        delete speciestree;

        //if (ctree)
        //    delete ctree;


    }



    //qDebug()<<QString::fromStdString(outstr);
    qDebug()<<qoutstr;

    qDebug()<<nbNull<<" null trees";

    //return a.exec();
}*/


