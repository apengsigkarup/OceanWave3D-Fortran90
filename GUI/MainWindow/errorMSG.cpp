#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMessageBox>
void MainWindow::error_grid_checkDomainLength()
{
    QMessageBox msgBox;
    msgBox.setText("Please check that the domain (length, width and depth) is consistent with grid file.\nNote; depth is the depth at the left boundary.");
    msgBox.exec();
}

void MainWindow::error_matlab()
{
    QMessageBox msgBox;
    msgBox.setText("Sorry, this version is compiled without support for MATLAB. \nFor MATLAB support please contact bo.paulsen@deltares.nl.\n\nNote; MATLAB scripts for post-processing are distributed with the source code.");
    msgBox.exec();
}


void MainWindow::errorMsgParameterCheck(){

    QMessageBox msgBox;
    msgBox.setText("Error checks are only implemented for:\n\nStream function\nJONSWAP\nP-M");
    msgBox.exec();

}

void MainWindow::errorMsgSFhZero(){

    QMessageBox msgBox;
    msgBox.setText("Incorrect water depth for stream function wave; h=0");
    msgBox.exec();

}

void MainWindow::errorMsgUnknownFile(){

    QMessageBox msgBox;
    msgBox.setText("Unknown input file");
    msgBox.exec();


}
