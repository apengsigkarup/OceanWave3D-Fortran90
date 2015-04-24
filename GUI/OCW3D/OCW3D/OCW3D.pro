#-------------------------------------------------
#
# Project created by QtCreator 2015-03-22T19:30:47
#
#-------------------------------------------------

QT       += core gui


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = OCW3D
TEMPLATE = app

SOURCES += main.cpp\
           mainwindow.cpp \
           convert.cpp

HEADERS  += mainwindow.h \
            convert.h \


FORMS    += mainwindow.ui

CONFIG += link_pkgconfig
PKGCONFIG += x11

INCLUDEPATH += /home/paulsen/MATLAB2014b/extern/include/
LIBS += -L /home/paulsen/MATLAB2014b/bin/glnxa64/
LIBS += -lmat
LIBS += -lmx



#LD_LIBRARY_PATH=/home/paulsen/MATLAB2014b/bin/glnxa64:/home/paulsen/MATLAB2014b/sys/os/glnxa64:$LD_LIBRARY_PATH

