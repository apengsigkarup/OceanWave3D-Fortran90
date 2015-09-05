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
           convert.cpp \
          ../QTwidgets/qcustomplot.cpp \
          customgrid.cpp \
          MainWindow/gridFunctions.cpp \
          MainWindow/aboutFunctions.cpp \
          MainWindow/runAndWrite.cpp \
          MainWindow/postProcessing.cpp \
          MainWindow/waveGeneration.cpp \
          MainWindow/errorMSG.cpp \
          MainWindow/check.cpp \
          checkdialog.cpp \


HEADERS  += mainwindow.h \
            convert.h \
            ../QTwidgets/qcustomplot.h \
            customgrid.h \
            checkdialog.h \
            versions.h \
            versions.h \




FORMS    += mainwindow.ui \
            checkdialog.ui \



CONFIG += link_pkgconfig console

#if _WIN32

#else
PKGCONFIG += x11
#endif



#if _WIN32

#INCLUDEPATH += ../QTwidgets/ \
#               C:\Program Files\boost\boost\
#else
INCLUDEPATH += /home/paulsen/MATLAB2014b/extern/include/ \
               ../QTwidgets/ \



#LIBS += -L/home/paulsen/lib/
#LIBS += -lOceanWave3D.a
#endif

#LIBS += -L/home/paulsen/MATLAB2014b/bin/glnxa64/
#LIBS += -lmat
#LIBS += -lmx




#LD_LIBRARY_PATH=/home/paulsen/MATLAB2014b/bin/glnxa64:/home/paulsen/MATLAB2014b/sys/os/glnxa64:$LD_LIBRARY_PATH

RESOURCES += \
    resources.qrc




#INCLUDEPATH += $$HOME/lib/
#DEPENDPATH += $$HOME/lib/

#unix:!macx:!symbian: PRE_TARGETDEPS += $$HOMElib/libOceanWave3D.a

