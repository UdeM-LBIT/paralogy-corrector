#############################################################################
# Makefile for building: ParalogyCorrector
# Generated by qmake (1.07a) (Qt 3.3.8b) on: Mon Jan  5 13:18:24 2015
# Project:  ParalogyCorrector.pro
# Template: app
# Command: $(QMAKE) "CONFIG+=RELEASE" -o Makefile ParalogyCorrector.pro
#############################################################################

####### Compiler, tools and options

CC       = gcc
CXX      = g++
LEX      = flex
YACC     = yacc
CFLAGS   = -pipe -Wall -W -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic  -D"NO_ENSEMBL=1" -DQT_NO_DEBUG -DQT_SHARED -DQT_TABLET_SUPPORT -DQT_THREAD_SUPPORT
CXXFLAGS = -pipe -std=c++0x -Wall -W -O2 -g -pipe -Wall -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector --param=ssp-buffer-size=4 -m64 -mtune=generic  -D"NO_ENSEMBL=1" -DQT_NO_DEBUG -DQT_SHARED -DQT_TABLET_SUPPORT -DQT_THREAD_SUPPORT
LEXFLAGS = 
YACCFLAGS= -d
INCPATH  = -I/usr/lib64/qt-3.3/mkspecs/default -I. -Isrc -I$(QTDIR)/include
LINK     = g++
LFLAGS   = 
LIBS     = $(SUBLIBS) -L$(QTDIR)/lib -ldl -lqt-mt -lXext -lX11 -lm
AR       = ar cqs
RANLIB   = 
MOC      = $(QTDIR)/bin/moc
UIC      = $(QTDIR)/bin/uic
QMAKE    = qmake
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY)
COPY_DIR = $(COPY) -r
INSTALL_FILE= $(COPY_FILE)
INSTALL_DIR = $(COPY_DIR)
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rmdir
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p

####### Output directory

OBJECTS_DIR = ./

####### Files

HEADERS = src/trees/treeiterator.h \
		src/trees/treeinfo.h \
		src/trees/polysolver.h \
		src/trees/node.h \
		src/trees/newicklex.h \
		src/div/util.h \
		src/div/define.h \
		src/trees/paralogycorrector.h
SOURCES = main.cpp \
		src/trees/treeiterator.cpp \
		src/trees/treeinfo.cpp \
		src/trees/polysolver.cpp \
		src/trees/node.cpp \
		src/trees/newicklex.cpp \
		src/trees/paralogycorrector.cpp
OBJECTS = main.o \
		treeiterator.o \
		treeinfo.o \
		polysolver.o \
		node.o \
		newicklex.o \
		paralogycorrector.o
FORMS = 
UICDECLS = 
UICIMPLS = 
SRCMOC   = 
OBJMOC = 
DIST	   = ParalogyCorrector.pro
QMAKE_TARGET = ParalogyCorrector
DESTDIR  = 
TARGET   = ParalogyCorrector

first: all
####### Implicit rules

.SUFFIXES: .c .o .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(UICDECLS) $(OBJECTS) $(OBJMOC)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJMOC) $(OBJCOMP) $(LIBS)

mocables: $(SRCMOC)
uicables: $(UICDECLS) $(UICIMPLS)

$(MOC): 
	( cd $(QTDIR)/src/moc && $(MAKE) )

Makefile: ParalogyCorrector.pro  /usr/lib64/qt-3.3/mkspecs/default/qmake.conf /usr/lib64/qt-3.3/lib/libqt-mt.prl
	$(QMAKE) "CONFIG+=RELEASE" -o Makefile ParalogyCorrector.pro
qmake: 
	@$(QMAKE) "CONFIG+=RELEASE" -o Makefile ParalogyCorrector.pro

dist: 
	@mkdir -p .tmp/ParalogyCorrector && $(COPY_FILE) --parents $(SOURCES) $(HEADERS) $(FORMS) $(DIST) .tmp/ParalogyCorrector/ && ( cd `dirname .tmp/ParalogyCorrector` && $(TAR) ParalogyCorrector.tar ParalogyCorrector && $(GZIP) ParalogyCorrector.tar ) && $(MOVE) `dirname .tmp/ParalogyCorrector`/ParalogyCorrector.tar.gz . && $(DEL_FILE) -r .tmp/ParalogyCorrector

mocclean:

uiclean:

yaccclean:
lexclean:
clean:
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) $(TARGET)


FORCE:

####### Compile

main.o: main.cpp src/trees/paralogycorrector.h \
		src/trees/newicklex.h \
		src/div/util.h \
		src/trees/node.h \
		src/div/define.h \
		src/trees/treeiterator.h \
		src/trees/treeinfo.h

treeiterator.o: src/trees/treeiterator.cpp src/trees/treeiterator.h \
		src/trees/node.h \
		src/div/define.h \
		src/trees/treeinfo.h \
		src/div/util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o treeiterator.o src/trees/treeiterator.cpp

treeinfo.o: src/trees/treeinfo.cpp src/trees/treeinfo.h \
		src/trees/node.h \
		src/div/util.h \
		src/div/define.h \
		src/trees/treeiterator.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o treeinfo.o src/trees/treeinfo.cpp

polysolver.o: src/trees/polysolver.cpp src/trees/polysolver.h \
		src/trees/node.h \
		src/trees/treeiterator.h \
		src/div/define.h \
		src/trees/treeinfo.h \
		src/div/util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o polysolver.o src/trees/polysolver.cpp

node.o: src/trees/node.cpp src/trees/node.h \
		src/div/define.h \
		src/trees/treeiterator.h \
		src/trees/treeinfo.h \
		src/div/util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o node.o src/trees/node.cpp

newicklex.o: src/trees/newicklex.cpp src/trees/newicklex.h \
		src/trees/node.h \
		src/div/util.h \
		src/div/define.h \
		src/trees/treeiterator.h \
		src/trees/treeinfo.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o newicklex.o src/trees/newicklex.cpp

paralogycorrector.o: src/trees/paralogycorrector.cpp src/trees/paralogycorrector.h \
		src/trees/node.h \
		src/trees/newicklex.h \
		src/div/define.h \
		src/trees/treeiterator.h \
		src/trees/treeinfo.h \
		src/div/util.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o paralogycorrector.o src/trees/paralogycorrector.cpp

####### Install

install:  

uninstall:  
