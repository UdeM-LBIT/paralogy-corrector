ParalogyCorrector
==============

Given a rooted binary gene tree + species tree, and a list of orthologous gene pairs, the program 
finds a gene tree that has the given orthology relationships (according to LCA reconciliation) 
that minimizes the Robinson-Foulds distance with the original tree.

The input is given as a .parc file, which can contain a list of trees to correct.
The format of the .parc file is a pseudo-xml format (see test.parc for an example).
The file has no root node and tag opening/closing must be preceded and followed by a newline character..
Each tree to correct is contained in a
```
<INSTANCE>
...
</INSTANCE>
tag.  Each instance tag must contain the child tags : 
<TREEID>
[The id of the gene tree.  Will be included in the output.]
</TREEID>

<GENETREE>
[The gene tree in newick format.  Each leaf must have a label, which will be used for the gene species mapping.]
</GENETREE>

<SPECIESTREE>
[The species tree in newick format.  Each leaf must have a label, which will be used for the gene species mapping.]
</SPECIESTREE>

<GENESPECIESMAPPING>
[
The map from the genes (leaf labels of the gene tree) to their corresponding species (leaf labels of the species tree).  
Each gene must be mapped.  
One mapping per line is given, and the format is
gene_name:species_name
]
</GENESPECIESMAPPING>

<CONSTRAINTS>
[
The gene orthologies to satisfy, using the leaf labels from the gene tree.
One constraint per line is given, and the format is
gene_name1:gene_name2
The given constraints do not need to be symmetric (ie gene_name2:gene_name1 does not need to be included - the output is unmodified whether it is included or not).
]
</CONSTRAINTS>
```

The output is a pseudo-xml with the following fields (one instance per tree): 
<INSTANCE>
<TREEID>
[the tree id of the current instance, as given in the .parc file]
</TREEID>
<BEFORE>
[the newick of the original gene tree, as given in the .parc file]
</BEFORE>
<AFTER>
[
the newick of the corrected gene tree, as given in the .parc file]
</AFTER>
[if an error occurred, the tag
<ERROR>
message
</ERROR>
will be added.  It typically comes up when all given orthologs were already orthologs in the given gene tree, and nothing was corrected.
]
</INSTANCE>



Example calls :

```
ParalogyCorrector -i "test.parc" -o "test.parc.out"
```


Command line arguments
----------------------

-i [infile] : infile is the .parc file that contains the required trees to correct/orthologies

-o [outfile]: (optional) the output is written in outfile, or on stdout if no outfile is specified

-f [old|new] : (optional, default=new) for backwards compatibility, setting -f old formats the output as in the older versions
