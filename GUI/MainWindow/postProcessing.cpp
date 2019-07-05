#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "advForce.h"

// Make a system call to gnuplot and print *.fort
//
void MainWindow::gnuplot(){

    QString fileName = ui->gnuplotFile->text();
    QString argument("/usr/bin/gnuplot -persist -e \"set grid;set autoscale fix;set xlabel \\\"x [m]\\\";set ylabel \\\" Free surface elevation, eta [m]\\\";plot \\\"" + fileName + "\\\" u 1:3 every ::2 w l;\" ");
    std::system(qPrintable(argument));

}

// Open a file dialog and slect a file *.fort to plot using gnuplot
//
void MainWindow::openFileDialog(){

    QString fileName = QFileDialog::getOpenFileName(
                this,
                "Select file",
                dir.currentPath(),
                "OCW3D (fort.*);;"
                );
    ui->gnuplotFile->setText(fileName);
    gnuplot();

}

// Read Kinematics.bin files
//
void MainWindow::readKinematicFile(){
    convertFiles.clearMem();
    ui->morisonX0->clear();
    ui->etaX0->clear();
    ui->convertStatus->clear();ui->convertStatus->repaint();
    ui->SelectOutput->setCurrentIndex(0);
    ui->readProgressBar->setVisible(true);
    convertFiles.read(ui->selectedPPfiles->text(),ui->readProgressBar);
    ui->readProgressBar->setValue(100);
    ui->convert->setEnabled(true);
    ui->SelectOutput->setEnabled(true);


}


//
//
void MainWindow::convertTo_setup(int output){

    if (output==1){ // morison force
        for (int i=0;(unsigned)i<convertFiles.x.size();i++){
            ui->morisonX0->addItem(QString::number(convertFiles.x[i]));
        }
        ui->morison_widget->setVisible(true);
        ui->eta_widget->setVisible(false);
    } else if (output==3){
        ui->etaX0->addItem("All locations");
        for (int i=0;(unsigned)i<convertFiles.x.size();i++){
            ui->etaX0->addItem(QString::number(convertFiles.x[i]));
        }


        ui->eta_widget->setVisible(true);
        ui->morison_widget->setVisible(false);


    } else {
        ui->morison_widget->setVisible(false);
        ui->eta_widget->setVisible(false);
    }


}

//
//
void MainWindow::convertTo(){
    ui->convertStatus->setVisible(true);
    ui->convertStatus->setText("Converting...");
    ui->convertStatus->repaint();

    int output = ui->SelectOutput->currentIndex();

    if (output==1){ // save morison Force
        ui->convertStatus->setText("Computing Morison forces");ui->convertStatus->repaint();
        convertFiles.force(ui->morisonX0->currentIndex(),ui->morison_D->value(),ui->Density->value(),ui->morison_cd->value(),ui->morison_cm->value());
        ui->convertStatus->setText("Done; Morison forces saved to working directory");
    }
    if (output==2){
#if MATLAB>0
        convertFiles.matlab();
#else
        error_matlab();
#endif
    }
    if (output==3){
        convertFiles.ascii(ui->selectedPPfiles->text(),ui->etaX0->currentIndex()-1);
        ui->convertStatus->setText("Done; Ascii file saved to working directory");
    }

#if externalOutputClass
    if (output==4){

        externalOutput *newOutputClass = new externalOutput(convertFiles,this);

        newOutputClass->show();
         ui->convertStatus->setText("Done; Input file to external code saved in working directory");

    }
#endif

}

void MainWindow::advancedMorison(){
    advForce*newOutputClass = new advForce(convertFiles,this);

    newOutputClass->show();


}


void MainWindow::selectPPfile(){

    QString fileName = QFileDialog::getOpenFileName(
                this,
                "Select file",
                dir.currentPath(),
                "OCW3D (Kinematics*.bin);;"
                );
    // We change working directory to the directory of "fileName"
    //
    dir.setCurrent(QFileInfo(fileName).path());
    ui->selectedPPfiles->setText(fileName);
    ui->workingDir->setText(QFileInfo(fileName).path());
}

//
//
