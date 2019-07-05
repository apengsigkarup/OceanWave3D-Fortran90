#include "mainwindow.h"
#include "ui_mainwindow.h"


void MainWindow::WaveTypeSelected(){

    int state = ui->DropDownListOutputType->currentIndex();

    if (state==0){
        ui->outputTable->setEnabled(false);
    }
    else {
        ui->outputTable->setEnabled(true);
        on_outputWidgetChanged(ui->nOutFiles->value());
    }

}


void MainWindow::storeASCII(bool checked){
    if (checked){
        ui->ACCII_label->setVisible(true);ui->ASCII_label2->setVisible(true);ui->nASCII->setVisible(true);}
    else {ui->ACCII_label->setVisible(false);ui->ASCII_label2->setVisible(false);ui->nASCII->setVisible(false);
    }
}

//
//
void MainWindow::on_outputWidgetChanged(int nFiles)
{
    int numberOfCol;

    // number of colums needed. For kinematics we need 6, for wave gagues we need 2
//    if ((ui->DropDownListOutputType->currentIndex()==1)&(ui->tableWidget->columnCount()!=6)) {
          if ((ui->DropDownListOutputType->currentIndex()==1)) {
        numberOfCol = 6;
        ui->tableWidget->setColumnCount(numberOfCol);ui->tableWidget->setRowCount(0);
        ui->tableWidget->setHorizontalHeaderLabels(QString("").split(""));
        ui->tableWidget->setHorizontalHeaderLabels(QString("Xmin;Xmax;Ymin;Ymax;tmin;tmax").split(";"));
        ui->tableWidget->setFixedSize(625,227);
    } else if ((ui->DropDownListOutputType->currentIndex()==2)) {
        numberOfCol = 2;
        ui->tableWidget->setColumnCount(numberOfCol);ui->tableWidget->setRowCount(0);
        ui->tableWidget->setHorizontalHeaderLabels(QString("x;y").split(";"));
        ui->tableWidget->setFixedSize(225,227);
    }



    QDoubleSpinBox *sp;
    if ((nFiles>0)&(~ui->tableWidget->isVisible())){
        ui->tableWidget->setVisible(true);
    }

    if (nFiles==0){
        ui->tableWidget->setVisible(false);
    }

    if (nFiles>ui->tableWidget->rowCount()){
        for (int i =ui->tableWidget->rowCount();i<nFiles;i++){

            ui->tableWidget->insertRow(i);

            for (int j = 0; j < ui->tableWidget->columnCount() ; j++) {
                ui->tableWidget->setCellWidget(i,j,new QDoubleSpinBox(ui->tableWidget));
                sp = (QDoubleSpinBox*)ui->tableWidget->cellWidget(i,j);
                sp->setMaximum(99999999);
                if ((j==ui->tableWidget->columnCount()-2)&(j>2)){sp->setValue(tStart);}
                if ((j==ui->tableWidget->columnCount()-1)&(j>2)){sp->setValue(ui->timeDuration->value());}
            }

        }
    }

    if (nFiles<ui->tableWidget->rowCount()){
        for (int i=ui->tableWidget->rowCount();i>nFiles;i--){

            ui->tableWidget->removeRow(i-1);
        }

    }

}

QSize MainWindow::QTableWidgetSize(QTableWidget *t) {
    int w = t->verticalHeader()->width() + 23 ; // +4 seems to be needed
    for (int i = 0; i < t->columnCount(); i++)
        w += t->columnWidth(i); // seems to include gridline (on my machine)

    int h = t->horizontalHeader()->height();
    for (int i = 0; i < t->rowCount(); i++)
        h += t->rowHeight(i) + 4;

    return QSize(w, h);
}
