#ifndef CUSTOMGRID_H
#define CUSTOMGRID_H

#include <QTableWidget>
#include "versions.h"
#include <QVector>
#include <QWidget>
#include <QHBoxLayout>
#include <qcustomplot.h>
#include <QMainWindow>


class customGrid
{
public:
    customGrid();
    ~customGrid();

    void readTabel(QTableWidget *,int,int,bool);
    void showGrid(double h, double L, int nx, int nz,bool);
    void showGrid();
    void laplaceSmoothing();
    void writeGridToFile(QString);
    void generateGrid();
    double L;
    double d;
    QVector<double> h;



    // test
private:
    // private member functions

    void computeGradients(QVector<double>, QVector<double>);
    void computeZ();
    double depth(double);


    // private data

    QVector<double> xData;
    QVector<double> zData;
    QVector<double> x;
    QVector<double> z;

    QVector<double> hx;
    QVector<double> hxx;

    int nx;
    int nz;
    bool sz;



};

#endif // CUSTOMGRID_H
