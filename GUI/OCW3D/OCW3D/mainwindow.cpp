#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iostream>
#include <fstream>
#include <QDebug>
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QVector>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QAbstractTableModel>
#include <math.h>
#include <QProcess>
#include "convert.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
/////////////////////////////


    ui->waveType->addItem("Wave theory");
    ui->waveType->addItem("Stream function");
    ui->waveType->addItem("JONSWAP");
    ui->waveType->addItem("P-M");
    ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);ui->irrFileName->setVisible(false);ui->irrFileName_label->setVisible(false);

    QRect SFwidget_geo = ui->widget_SF->geometry();
    ui->widget_JONSWAP->setGeometry(SFwidget_geo);
    ui->widget_JONSWAP->setVisible(false);
    ui->aboutText_OCW3dGUI->setVisible(false);
    ui->aboutText_OCW3D_publications->setVisible(false);
    ui->readProgressBar->setVisible(false);

    ui->tableWidget->setColumnCount(6);
    ui->tableWidget->setHorizontalHeaderLabels(QString("Xmin;Xmax;Ymin;Ymax;tmin;tmax").split(";"));
    ui->tableWidget->setVisible(false);
    ui->ACCII_label->setVisible(false);ui->ASCII_label2->setVisible(false);ui->nASCII->setVisible(false);

    ui->SelectOutput->setEnabled(false);
    ui->convert->setEnabled(false);
    ui->morison_widget->setVisible(false);
    ui->convertStatus->setVisible(false);
    ui->convertWarning->setVisible(false);ui->convertWarning->viewport()->setAutoFillBackground(false);

    connect(ui->waveType,SIGNAL(currentIndexChanged(int)),this,SLOT(on_waveTheoryChanged(int)));
    connect(ui->storeAscii_onOff,SIGNAL(clicked(bool)),this,SLOT(storeASCII(bool)));
    connect(ui->LorP_ComboBox,SIGNAL(currentIndexChanged(int)),this,SLOT(on_waveTheoryChanged()));
    connect(ui->nOutFiles,SIGNAL(valueChanged(int)),this,SLOT(on_outputWidgetChanged(int)));
    connect(ui->OpenDirBrowser,SIGNAL(clicked()),this,SLOT(openWorkDirDialog()));
    connect(ui->selectPPfile,SIGNAL(clicked()),this,SLOT(selectPPfile()));
    connect(ui->run,SIGNAL(clicked()),this,SLOT(run()));
    connect(ui->selectGnuplotFile,SIGNAL(clicked()),this,SLOT(openFileDialog()));
    connect(ui->plot,SIGNAL(clicked()),this,SLOT(gnuplot()));
    connect(ui->read_bottom,SIGNAL(clicked()),this,SLOT(readKinematicFile()));
    connect(ui->about_combobox,SIGNAL(currentIndexChanged(int)),this,SLOT(about_changed(int)));
    connect(ui->SelectOutput,SIGNAL(currentIndexChanged(int)),this,SLOT(convertTo_setup(int)));
    connect(ui->convert,SIGNAL(clicked()),this,SLOT(convertTo()));
    ui->workingDir->setText(dir.currentPath());
//    qApp->setStyleSheet(QString("QWidget {background:rgb(0,0,250); olor gray;}"));
}

MainWindow::~MainWindow()
{
    delete ui;
}



void MainWindow::about_changed(int i){

    QRect pos = ui->aboutText_OCW3D->geometry();

    if (i==0){
        ui->aboutText_OCW3D->setVisible(true);
        ui->aboutText_OCW3D_publications->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(false);

    }
    if(i==1){
        ui->aboutText_OCW3dGUI->setGeometry(pos);
        ui->aboutText_OCW3D->setVisible(false);
        ui->aboutText_OCW3D_publications->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(true);
    }

    if (i ==2){
        ui->aboutText_OCW3D_publications->setGeometry(pos);
        ui->aboutText_OCW3D->setVisible(false);
        ui->aboutText_OCW3dGUI->setVisible(false);
        ui->aboutText_OCW3D_publications->setVisible(true);

    }



}


void MainWindow::run(){
    QProcess *OCW = new QProcess();
    OCW->start("/usr/bin/xterm -hold -e \\\"OceanWave3D\\\"");

//    system("/usr/bin/xterm -e /home/paulsen/bin/OceanWave3D");

}

void MainWindow::gnuplot(){
    QString fileName = ui->gnuplotFile->text();
    QString argument("/usr/bin/gnuplot -persist -e \"set grid;set autoscale fix;set xlabel \\\"x [m]\\\";set ylabel \\\" Free surface elevation, eta [m]\\\";plot \\\"" + fileName + "\\\" u 1:3 w l;\" ");
    std::system(qPrintable(argument));

}

void MainWindow::openFileDialog(){


    QString selfilter = tr("OCW3D (fort.*)");
    QString fileName = QFileDialog::getOpenFileName(
            this,
                "Select file",
            dir.currentPath(),
            tr("OCW3D (fort.*);" ),
            &selfilter
    );
    ui->gnuplotFile->setText(fileName);
    gnuplot();

}

void MainWindow::openWorkDirDialog(){

    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly);
    dialog.exec();

    dir.setCurrent(dialog.directory().path());

    ui->workingDir->setText(dialog.directory().path());
}

void MainWindow::readKinematicFile(){
    ui->convertWarning->setVisible(false);
    ui->readProgressBar->setVisible(true);
    convertFiles.read(ui->selectedPPfiles->text(),ui->readProgressBar);
    ui->readProgressBar->setValue(100);
    ui->convert->setEnabled(true);
    ui->SelectOutput->setEnabled(true);
}

void MainWindow::saveAsAscii(){

convertFiles.ascii();

}

void MainWindow::convertTo_setup(int output){

    if (output==1){ // morison force
        double dx = convertFiles.x[1]-convertFiles.x[0];
        ui->morison_x0->setMaximum(convertFiles.xend*dx);
        ui->morison_x0->setMinimum(convertFiles.xbeg*dx);
        ui->morison_widget->setVisible(true);
    } else {
        ui->morison_widget->setVisible(false);
    }

}
void MainWindow::convertTo(){
    ui->convertStatus->setVisible(true);
    ui->convertStatus->setText("Converting...");
    qApp->processEvents();
    int output = ui->SelectOutput->currentIndex();

    if (output==1){ // save morison Force
        convertFiles.force(ui->morison_x0->value(),ui->morison_D->value(),ui->Density->value(),ui->morison_cd->value(),ui->morison_cm->value());
    }
    if (output==2){
        convertFiles.matlab();
    }
    if (output==3){
        convertFiles.ascii();
    }
    ui->convertStatus->setText("DONE");
}
void MainWindow::saveAsMat(){

    convertFiles.matlab();
}
void MainWindow::selectPPfile(){
    PPfile.setDirectory(dir.currentPath());

    PPfile.exec();

//    std::cout << PPfile.selectedFiles().join("\n").toStdString() << std::endl;
    ui->selectedPPfiles->setText(PPfile.selectedFiles().join(";"));
    ui->convertWarning->setVisible(true);
}

void MainWindow::storeASCII(bool checked){
    if (checked){
        ui->ACCII_label->setVisible(true);ui->ASCII_label2->setVisible(true);ui->nASCII->setVisible(true);}
    else {ui->ACCII_label->setVisible(false);ui->ASCII_label2->setVisible(false);ui->nASCII->setVisible(false);
    }
}

void MainWindow::on_outputWidgetChanged(int nFiles)
{
     QSpinBox *sp;
    if (nFiles>0&~ui->tableWidget->isVisible()){
        ui->tableWidget->setVisible(true);
    }

    if (nFiles==0){
        ui->tableWidget->setVisible(false);
    }

    if (nFiles>ui->tableWidget->rowCount()){
        for (int i =ui->tableWidget->rowCount();i<nFiles;i++){

            std::cout << "add row " << i << std::endl;
            ui->tableWidget->insertRow(i);
            for (int j = 0; j < 6 ; j++) {
                   ui->tableWidget->setCellWidget(i,j,new QSpinBox(ui->tableWidget));
                  sp = (QSpinBox*)ui->tableWidget->cellWidget(i,j);
                  sp->setMaximum(99999999);
            }


        }
    }

    if (nFiles<ui->tableWidget->rowCount()){
        for (int i=ui->tableWidget->rowCount();i>nFiles;i--){
            std::cout << "removing row " << i << std::endl;
            ui->tableWidget->removeRow(i-1);
        }

    }

}

void MainWindow::on_waveTheoryChanged()
{
    MainWindow::on_waveTheoryChanged(1);
}

void MainWindow::on_waveTheoryChanged(int waveTheory)
{
    if (waveTheory==1) {
        ui->widget_SF->setVisible(true);ui->LorP_ComboBox->setVisible(true);

        if (ui->LorP_ComboBox->currentIndex()==0){
            ui->SF_TL_unit->setText("m");
            ui->SF_TLabel->setText("Wave length");
        }
        if (ui->LorP_ComboBox->currentIndex()==1) {
            ui->SF_TL_unit->setText("s");
            ui->SF_TLabel->setText("Wave period");
        }
        ui->widget_JONSWAP->setVisible(false);
    } else if (waveTheory==2){
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
        ui->widget_JONSWAP->setVisible(true);
    } else if (waveTheory==3){
        ui->widget_JONSWAP->setVisible(true);
        ui->widget_SF->setVisible(false);
    } else {
        ui->widget_JONSWAP->setVisible(false);
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
    }
}

void MainWindow::on_pushButton_2_clicked()
{
// Read inputs from tab GENERAL
    // Read header
    QString header = ui->header_input->toPlainText();
    // Read mode
    int mode;
    bool nLin = ui->Nonlin_onOff->isChecked();
    bool Lin = ui->linear_onOff->isChecked();
(nLin>Lin) ? mode = 1 : mode=0;
    // Gravity
    double gravity = ui->gravity_input->value();
    double Density = ui->Density->value();
    // Breaking filter
    double breaking = ui->breaking_beta0->value();
    if (breaking==0){breaking = 1000;}
    // Time duration
    double timeDur = ui->timeDuration->value();
    double tStart = ui->tStart->value();
    // Read initial conditions
    int initialConditions = ui->initialCondition->currentIndex();

   // read discretization
    int alpha = ui->alpha->value();
    int beta = ui->beta->value();
    int gamma = ui->gamma->value();
    int a = ui->a->value();
    int b = ui->b->value();
    int c = ui->c->value();
    double dt = ui->dt->value();

    // Compute the number of time steps
    int Nsteps = timeDur/dt;


    // Read geometry
    int length = ui->length ->value();
    int width = ui->width->value();
    int depth = ui->depth->value();
    int nx = ui->nx->value();
    int ny = ui->ny->value();
    int nz = ui->nz->value();
    bool sx = ui->sx->isChecked();
    bool sy = ui->sy->isChecked();
    bool sz = ui->sz->isChecked();

// Read wave generation
    // wave theory
    int waveTheory = ui->waveType->currentIndex();
    irrFilename = ui->irrFileName->text();
    if (waveTheory==1){ // Stream Fucntion
        SF_H = ui->SF_H->value();
        SF_T = ui->SF_T->value();
        SF_h = ui->SF_h->value();
        SF_U = ui->SF_U->value();
        SorE = ui->stokesOrEuler->currentIndex();
        SF_n = ui->SF_n->value();
        LorT = ui->LorP_ComboBox->currentIndex();
        SF_L=100;//Dummy value
       if (LorT==0) {SF_L =SF_T;}
    }
    if (waveTheory==2||waveTheory==3){ //JONSWAP or PM
        Hs = ui->Hs->value();
        Tp = ui->Tp->value();
        h = ui->h->value();
        gamma = ui->gamma->value();
        seed = ui->seed->value();
        khmax = ui->maxkh->value();
        i_spec;
        (waveTheory==2) ? i_spec=1 : i_spec=0;
    }

    // Read relaxation / damping zones
        //Generation
        double xGen[2];double yGen[2];
        xGen[0] = ui->xGenStart->value();xGen[1] = ui->xGenEnd->value();
        yGen[0] = ui->yGenStart->value();yGen[1] = ui->yGenEnd->value();
        double rampTime = ui->rampTime->value();

        double xAbsorb[2];double yAbsorb[2];
        xAbsorb[0] = ui->xAbsorbStart->value();xAbsorb[1] = ui->xAbsorbEnd->value();
        yAbsorb[0] = ui->yAbsorbStart->value();yAbsorb[1] = ui->yAbsorbEnd->value();


    // Read output files
    bool storeAscii = ui->storeAscii_onOff->isChecked();
    int nASCIIsteps = ui->nASCII->value();

    // compute resolution
    double epsilon = 1e-12;
    double dx = length/(nx-1);
    double dy = (width+epsilon)/(ny-1);
    double resolution[] ={dx,dx,dy,dy,dt,dt};
    QSpinBox* sp;
    double outputValues[ui->tableWidget->rowCount()][ui->tableWidget->columnCount()];
    for (int i=0;i<ui->tableWidget->rowCount();i++){
        for (int j = 0; j < 6 ; j++) {
            sp = (QSpinBox*)ui->tableWidget->cellWidget(i,j);
            outputValues[i][j]= round(sp->value()/resolution[j]);

        }
        if (ny==1) {outputValues[i][2]=1;outputValues[i][3]=1;}
        if (outputValues[i][5]>Nsteps){outputValues[i][5]=Nsteps;}
        if (outputValues[i][4]==0){outputValues[i][4]=1;}

     }

//     Write OceanWave3D input file
    QFile myfile("OceanWave3D.inp");
    myfile.open (QIODevice::WriteOnly);
    QTextStream outStream(&myfile);
    outStream << header << "\n";
    outStream << initialConditions << " " << waveTheory << " " << breaking << "\n";
    outStream << length << " "<< width << " " << depth << " " << nx << " " << ny << " " << nz << " " << sx << " " << sy << " " << sz << " 1 1 1\n";
    outStream << alpha << " " << beta << " " << gamma << " " << a << " " << b << " " << c << "\n";
    outStream << Nsteps << " " << dt << " " << "1 0 0 " << tStart << "\n";
    outStream << gravity <<" " << Density << "\n";
    outStream << "1 1 0 23 1e-8 1e-6 1 V 1 1 2 \n";
    if (waveTheory==1){
        outStream << SF_H << " " << SF_h << " " << SF_L << " " << SF_T << " "<< LorT << " " << SF_U << " " << SorE << " " << "8 " << " " << SF_n << "\n";
    } else { outStream << "1 1 1 1 1 1 1 1 1\n";}

    if (ui->nOutFiles->value()==0&~storeAscii){ // No outputs requested
        outStream << "0 0 \n";
    } else { // Outputs are requested

        if (ui->nOutFiles->value()>0&storeAscii){ // Both binary and ascii files
            outStream << -1*nASCIIsteps << " 20 0 " << ui->nOutFiles->value() << "\n";
            }
        if (ui->nOutFiles->value()==0&storeAscii){ // Only ASCII files
            outStream << -1*nASCIIsteps << " 1 0 0\n";
            }
        if (ui->nOutFiles->value()>0&~storeAscii){ // Only binary Files
            outStream << "1 20 1 " << ui->nOutFiles->value() << "\n";
        }

        if (ui->nOutFiles->value()>0){ // Writing output positions for binary files

            for (int i = 0;i<ui->nOutFiles->value();i++){ // loop over rows
                for (int j = 0; j<6;j=j+2){
                    outStream << outputValues[i][j] << " " << outputValues[i][j+1] << " 1 ";
                }
                outStream << "\n";
            }

        }

   }

    outStream << mode << " 0\n"; // We hardcode the surface pressure term
    outStream << "0 6 10 0.08 0.08 0.4 \n"; // SG-filtering
    if (ui->pressureDampingOrRelax->currentIndex()==0) {
        outStream << "1 " << rampTime << " 1 X 0 \n";
        outStream << xGen[0] << " " << xGen[1] << " " << yGen[0] << " " << yGen[1] << " 10 3.5 X 1 X 0 \n";
        outStream << "1 1\n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 1 1 0\n";
    }
    if (ui->pressureDampingOrRelax->currentIndex()==1) {
        outStream << "1 " << rampTime << " 2 X 0 \n";
        outStream << xGen[0] << " " << xGen[1] << " " << yGen[0] << " " << yGen[1] << " 10 3.5 X 1 X 0 \n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 1 9 3.5 X 0 X 0 \n";
        outStream << "0 0\n";
    }

    outStream << "0 0 0 0 0 0 0\n"; // SWENSE
    outStream << "0 \n";
    // Irregualr waves
    if (waveTheory>1){outStream << i_spec << " " << Tp << " " << Hs << " " << depth << " " << khmax << " " << seed << " " << seed << " 0 0 " << irrFilename << " 3D-off1 0\n";}
    //
    myfile.close();
    ui->run->setEnabled(true);
}
