#ifndef advForce_H
#define advForce_H

#include <QWidget>
#include "convert.h"
#include <QGridLayout>
#include <QGroupBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QFont>
#include <QTableWidget>
#include "qcustomplot.h"
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QSpacerItem>
#include <QFile>
#include <QFileInfo>
#include <QFileDialog>
#include <QComboBox>


class advForce : public QWidget
{
    Q_OBJECT
public:
    explicit advForce(const convert& data_, QWidget *parent = 0);
    ~advForce();
    const convert *data;
    void interp_linear(const std::vector<double> f_ind,const std::vector<double> x_ind,const std::vector<double> x_out, std::vector<double>& f_out, int& id);
    void computeForce(int nn);
signals:
    
public slots:

private:

    // functions
    void createLeftSide();
    void createRightSide();
    void updateTableSize(QTableWidget *, int, int);
    void readGeometryFile();
    void readTable();
    void saveGeometry();

    QVector<double> Fout;
    std::vector<double> D_table,Cm_table,Cd_table,z_table;

    QFileInfo selectFile();

    QWidget *ui;

    QGridLayout *mainLayout;

    QGroupBox *leftView;
    QGroupBox *leftView2;
    QGroupBox *leftView3;
    QGroupBox *rightView;

    QLineEdit *advForceFile;
    QLineEdit *geometryFile;
    QLineEdit *geometryHeader;

    QPushButton *selectadvForceFile;
    QPushButton *selectGeometryFile;
    QPushButton *updatePlotAndFile;
    QPushButton *save;
    QPushButton *plot;


    QLabel *density_label;
    QLabel *selectadvForceFile_lable;
    QLabel *nradvForceLevels_label;
    QLabel *nVerticalLayers_label;
    QLabel *selectGeometryFile_label;
    QLabel *xLocation_label;
    QLabel *geometryHeader_label;

    QSpinBox *nradvForceLevels;
    QSpinBox  *nRows;

    QDoubleSpinBox *density;

    QComboBox *xLocation;
    QComboBox *gridType;

    QFileInfo advForceFile_fileInfo;
    QFileInfo geometryFile_fileInfo;

    QTableWidget *forceCoeffsAndDiameter;

    QCustomPlot *monopilePlot;

    std::vector<double> F;
    std::vector<double> F_advForce;


private slots:
    void updateTableSize_();
    void setadvForceFile_();
    void setGeometryFile_();
    void plotMonoPile_();
    void save_();
    void plot_();
    void gridTypeChanged_();

};

#endif // advForce_H
