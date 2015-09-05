#include "mainwindow.h"
#include "checkdialog.h"
#include <QApplication>
#include "fstream"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
//     QCoreApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
