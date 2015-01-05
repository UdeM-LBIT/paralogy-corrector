#-------------------------------------------------
#
# Project created by QtCreator 2013-05-15T13:59:22
#
#-------------------------------------------------

INCLUDEPATH += src
DEPENDPATH += src

QT       += core

QT       -= gui

TARGET = ParalogyCorrector
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp \
    src/trees/treeiterator.cpp \
    src/trees/treeinfo.cpp \
    src/trees/polysolver.cpp \
    src/trees/node.cpp \
    src/trees/newicklex.cpp \
    src/trees/paralogycorrector.cpp

HEADERS += \
    src/trees/treeiterator.h \
    src/trees/treeinfo.h \
    src/trees/polysolver.h \
    src/trees/node.h \
    src/trees/newicklex.h \
    src/div/util.h \
    src/div/define.h \
    src/trees/paralogycorrector.h




unix{

    LIBS += -ldl
}

win32{

}


#QT += sql
DEFINES += "NO_ENSEMBL=1"

QMAKE_CXXFLAGS += -std=c++0x
