#include "mainwindow.h"
#include "ui_mainwindow.h"

// Update about tab
void MainWindow::about_changed(int i){

    QRect pos = ui->aboutText_OCW3D->geometry();

    if (i==0){
        ui->aboutText_OCW3D->setVisible(true);
        ui->aboutText_OCW3D_publications->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(false);
        ui->aboutText_OCW3DVersion->setVisible(false);
    }
    if (i==1){
        ui->aboutText_OCW3dGUI->setGeometry(pos);
        ui->aboutText_OCW3D->setVisible(false);
        ui->aboutText_OCW3D_publications->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(true);
        ui->aboutText_OCW3DVersion->setVisible(false);
    }

    if (i==2){
        ui->aboutText_OCW3D_publications->setGeometry(pos);
        ui->aboutText_OCW3D->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(false);
        ui->aboutText_OCW3D_publications->setVisible(true);
        ui->aboutText_OCW3DVersion->setVisible(false);
    }
    if (i==3){
        ui->aboutText_OCW3DVersion->setGeometry(pos);
        ui->aboutText_OCW3D->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(false);
        ui->aboutText_OCW3D_publications->setVisible(false);
        ui->aboutText_OCW3DVersion->setVisible(true);
    }



}
