#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "convert.h"
#include "mat.h"
#include <QFileDialog>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButton_2_clicked();
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
    void saveAsAscii(); //remove
    void saveAsMat(); // remove
    void convertTo_setup(int);
    void convertTo();


private:
    Ui::MainWindow *ui;
    QFileDialog dialog;
    QFileDialog PPfile;
    QFileDialog gnuplotFile;

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
    double khmax;
    int SF_n;
    int LorT;
    int seed;
    int i_spec;
    QString irrFilename;
    QDir dir;
    convert convertFiles;
};

#endif // MAINWINDOW_H
