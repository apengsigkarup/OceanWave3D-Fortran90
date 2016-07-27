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
#include <QMessageBox>
#include <QAbstractTableModel>
#include <math.h>
#include <QProcess>

#include "convert.h"
extern "C" void MAIN_();

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    checkDialog c;
    c.show();

    // shortcuts
    ui->actionOpen->setShortcut(Qt::Key_O | Qt::CTRL);
    ui->actionChange_directory->setShortcut(Qt::Key_D | Qt::CTRL);
    ui->actionClose->setShortcut(Qt::Key_Q | Qt::CTRL);
    ui->actionRun_F5->setShortcut(Qt::Key_F5);
    ui->actionClear_case_F6->setShortcut(Qt::Key_F6);
    ui->actionWrite_input_file_F4->setShortcut(Qt::Key_F4);
    ui->actionCheck->setShortcut(Qt::Key_F7);
    QSignalMapper *signalMapper = new QSignalMapper(this);
    for( int index=0; index < ui->output->count() ; ++index ){
        QShortcut *shortcut  =new QShortcut( QKeySequence(QString("Ctrl+%1").arg( index +1 ) ), this );
        connect( shortcut, SIGNAL(activated() ), signalMapper, SLOT( map() ) );
        signalMapper->setMapping( shortcut, index );
    }
    connect( signalMapper, SIGNAL(mapped( int )),ui->output, SLOT(setCurrentIndex( int )) );


    // Wave theory
    ui->widget_SF->setVisible(false);
    ui->LorP_ComboBox->setVisible(false);
    ui->widget_waveFile->setVisible(false);
    QRect SFwidget_geo = ui->widget_SF->geometry();
    ui->widget_JONSWAP->setGeometry(SFwidget_geo);
    ui->widget_JONSWAP->setVisible(false);
    ui->widget_waveFile->setGeometry(SFwidget_geo);
    ui->widget_customSpectrum->setVisible(false);
    ui->widget_customSpectrum->setGeometry(SFwidget_geo);

    // About
    ui->aboutText_OCW3dGUI->setVisible(false);
    ui->aboutText_OCW3D_publications->setVisible(false);
    ui->aboutText_OCW3DVersion->setVisible(false);

    // Post processing
    ui->readProgressBar->setVisible(false);
    ui->tableWidget->setColumnCount(6);
    ui->tableWidget->setVisible(false);
    ui->SelectOutput->setEnabled(false);
    ui->convert->setEnabled(false);
    ui->morison_widget->setVisible(false);
    ui->eta_widget->setVisible(false);
    QRect Morison_geo = ui->morison_widget->geometry();
    ui->eta_widget->setGeometry(Morison_geo);
    ui->convertStatus->setVisible(false);

    // Custom grid
    ui->geometry_table->setRowCount(2);
    ui->geometry_table->setVerticalHeaderLabels(QString("x [m];d [m]").split(";"));
    ui->customGridWidget->setVisible(false);

    // set labels
    ui->alpha_label->setText(QString((QChar) 0x03B1));
    ui->beta_label->setText(QString((QChar) 0x03B2));
    ui->gamma_label->setText(QString((QChar) 0x03B3));

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
    connect(ui->showGrid,SIGNAL(clicked()),this,SLOT(showGrid()));
    connect(ui->geometryType,SIGNAL(currentIndexChanged(int)),this,SLOT(geometryType_changed(int)));
    connect(ui->nGridPoints,SIGNAL(valueChanged(int)),this,SLOT(nGridPoints_changed()));
    connect(ui->smooth,SIGNAL(clicked()),this,SLOT(smooth()));
    connect(ui->generateGrid,SIGNAL(clicked()),this,SLOT(generateGrid()));
    connect(ui->selectWaveFile,SIGNAL(clicked()),this,SLOT(selectWaveFile()));
    connect(ui->selectGridFile,SIGNAL(clicked()),this,SLOT(selectGridFile()));
    connect(ui->selectWaveFile_eta,SIGNAL(clicked()),this,SLOT(selectWaveFile_eta()));
    connect(ui->DropDownListOutputType,SIGNAL(currentIndexChanged(int)),this,SLOT(WaveTypeSelected()));
    connect(ui->pushButton_advancedMorison,SIGNAL(clicked()),this,SLOT(advancedMorison()));
    // default
    connect(ui->checkBox_constantWidget,SIGNAL(stateChanged(int)),this,SLOT(constantWidget()));
    connect(ui->checkBox_breakingWidget,SIGNAL(stateChanged(int)),this,SLOT(breakingWidget()));
    connect(ui->checkBox_FD,SIGNAL(stateChanged(int)),this,SLOT(FDWidget()));
    connect(ui->checkBox_precon,SIGNAL(stateChanged(int)),this,SLOT(preconWidget()));

    // ations
    connect(ui->actionClose,SIGNAL(triggered()),this,SLOT(close()));
    connect(ui->actionOpen,SIGNAL(triggered()),this,SLOT(openFile()));
    connect(ui->actionRun_F5,SIGNAL(triggered()),this,SLOT(run()));
    connect(ui->actionClear_case_F6,SIGNAL(triggered()),this,SLOT(clearCase()));
    connect(ui->actionWrite_input_file_F4,SIGNAL(triggered()),this,SLOT(writeInputFile()));
    connect(ui->actionChange_directory,SIGNAL(triggered()),this,SLOT(openWorkDirDialog()));
    connect(ui->actionCheck,SIGNAL(triggered()),this,SLOT(checkCase()));
    ui->workingDir->setText(dir.currentPath());

    // Special versions

#if externalOutputClass
    ui->SelectOutput->addItem("External output");
#endif


}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::constantWidget(){

    if (ui->checkBox_constantWidget->isChecked()){
    ui->constantWidget->setEnabled(false);
    ui->Density->setValue(1025);
    ui->gravity_input->setValue(9.82);
    }else {ui->constantWidget->setEnabled(true);}
}
void MainWindow::breakingWidget(){

    if (ui->checkBox_breakingWidget->isChecked()){
    ui->breakingWidget->setEnabled(false);
    ui->breaking_beta0->setValue(0);
    }else {ui->breakingWidget->setEnabled(true);}
}

void MainWindow::FDWidget(){

    if (ui->checkBox_FD->isChecked()){
    ui->FDwidget->setEnabled(false);
    ui->alpha->setValue(3);
    ui->beta->setValue(3);
    ui->gamma->setValue(3);
    }else {ui->FDwidget->setEnabled(true);}
}
void MainWindow::preconWidget(){

    if (ui->checkBox_precon->isChecked()){
    ui->PreconWidget->setEnabled(false);
    ui->a->setValue(1);
    ui->b->setValue(1);
    ui->c->setValue(1);
    }else {ui->PreconWidget->setEnabled(true);}
}

#include "MainWindow/gridFunctions.cpp"
#include "MainWindow/aboutFunctions.cpp"
#include "MainWindow/runAndWrite.cpp"
#include "MainWindow/postProcessing.cpp"
#include "MainWindow/waveGeneration.cpp"


double MainWindow::dispersion_L(double h,double L){

    double k = 2*M_PI/L;
    return  sqrt(constants::g * k *tanh(k*h));

}

double MainWindow::dispersion_T(double h, double T){
    double omega = 2*M_PI/T;
    double x = h*omega/sqrt(constants::g*h);

    double kh = pow(x,2)*pow(1-exp(-pow(x,constants::beta)),-1/constants::beta);
    return kh/h;
}

double MainWindow::round(double x)
{
   return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
}





