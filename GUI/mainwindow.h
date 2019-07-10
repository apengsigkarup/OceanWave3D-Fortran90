#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "versions.h"
#include "math.h"
#include <cmath>
#include "checkdialog.h"
#include "convert.h"
#include <QFileDialog>
#include <QWidget>
#include <QSignalMapper>
#include <QShortcut>
#include <QHBoxLayout>
#include <QKeySequence>
#include "advForce.h"

#include "QTwidgets/qcustomplot.h"
#include "customgrid.h"


#if externalOutputClass
#include "MainWindow/externaloutput.h"
#endif

namespace Ui {
class MainWindow;
}

namespace constants
{
    const double g = 9.81;
    const double beta = 2.4908;
} // namespace constants

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void writeInputFile();
    void on_waveTheoryChanged(int waveTheory);
    void on_waveTheoryChanged();
    void on_outputWidgetChanged(int nFiles);
    void storeASCII(bool checked);
    void openFileDialog();
    void openWorkDirDialog();
    void selectPPfile();
    void run();
    void gnuplot();
    void readKinematicFile();
    void about_changed(int);
    void convertTo_setup(int);
    void convertTo();
    void showGrid();
    void smooth();
    void generateGrid();
    void geometryType_changed(int);
    void nGridPoints_changed();
    void openFile();
    void clearCase();
    void selectWaveFile();
    void selectWaveFile_eta();
    void selectGridFile();
    void WaveTypeSelected();
    void advancedMorison();


    // use constant
    void constantWidget();
    void breakingWidget();
    void FDWidget();
    void preconWidget();

    // error messages
    void error_grid_checkDomainLength();
    void error_matlab();
    void errorMsgParameterCheck();
    void errorMsgSFhZero();
    void errorMsgUnknownFile();
    // check case
    bool checkCase();


private:
    Ui::MainWindow *ui;
    QFileDialog dialog;
    double dispersion_T(double,double);
    double dispersion_L(double,double);
    double round(double);
    QSize QTableWidgetSize(QTableWidget*);
    void readAndSplit(QFile &,QStringList &);

    // variables
    double SF_H;
    double SF_T;
    double SF_h;
    double SF_U;
    double SF_L;
    double SorE;
    double Hs;
    double Tp;
    double h;
    double gamma;
    double gamma_jonswap;
    double khmax;
    int SF_n;
    int LorT;
    int seed;
    int i_spec;
    QString irrFilename;
    QDir dir;
    convert convertFiles;
    customGrid grid;
    double tStart;

};

#endif // MAINWINDOW_H
