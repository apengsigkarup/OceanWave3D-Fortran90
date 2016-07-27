#-------------------------------------------------
#
# Project created by QtCreator 2015-03-22T19:30:47
#
#-------------------------------------------------

QT       += core gui


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
TARGET = OCW3D
TEMPLATE = app

#    Release:DESTDIR = release
#    Release:OBJECTS_DIR = release/.obj
#    Release:MOC_DIR = release/.moc
#    Release:RCC_DIR = release/.rcc
#    Release:UI_DIR = release/.ui

#    Debug:DESTDIR = debug
#    Debug:OBJECTS_DIR = debug/.obj
#    Debug:MOC_DIR = debug/.moc
#    Debug:RCC_DIR = debug/.rcc
#    Debug:UI_DIR = debug/.ui


SOURCES += main.cpp\
           mainwindow.cpp \
           convert.cpp \
           QTwidgets/qcustomplot.cpp \
           customgrid.cpp \
           MainWindow/gridFunctions.cpp \
           MainWindow/aboutFunctions.cpp \
           MainWindow/runAndWrite.cpp \
           MainWindow/postProcessing.cpp \
           MainWindow/waveGeneration.cpp \
           MainWindow/errorMSG.cpp \
           MainWindow/check.cpp \
           checkdialog.cpp \
           MainWindow/externaloutput.cpp \
    MainWindow/outputLocations.cpp \
    advForce.cpp


HEADERS  += mainwindow.h \
            convert.h \
            QTwidgets/qcustomplot.h \
            customgrid.h \
            checkdialog.h \
            versions.h \
            versions.h \
    MainWindow/externaloutput.h \
advForce.h

FORMS    += mainwindow.ui \
            checkdialog.ui \
    MainWindow/externaloutput.ui

CONFIG += link_pkgconfig console debug_and_release


PKGCONFIG += x11




#if _WIN32

#INCLUDEPATH += ../QTwidgets/ \
#               C:\Program Files\boost\boost\
#else
INCLUDEPATH += /home/paulsen/MATLAB2014b/extern/include/ \
               QTwidgets/ \



#if MATLAB
#INCLUDEPATH += /home/paulsen/MATLAB2014b/extern/include/
#LIBS += -L/home/paulsen/MATLAB2014b/bin/glnxa64/
#LIBS += -lmat
#LIBS += -lmx
#end

RESOURCES += \
    resources.qrc
