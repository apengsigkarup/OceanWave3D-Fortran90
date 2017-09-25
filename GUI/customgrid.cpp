#include "customgrid.h"
#include <iostream>
#include <math.h>
#include <QDoubleSpinBox>

// Constructor
customGrid::customGrid()
{
}

// Destructor
customGrid::~customGrid(){
}

void customGrid::readTabel(QTableWidget *gridTable,int nx_,int nz_,bool sz_){

    QDoubleSpinBox *sp;
int n = gridTable->columnCount() ;
xData.resize(n);
zData.resize(n);

    for (int i = 0;i<n;i++){
        sp = (QDoubleSpinBox*)gridTable->cellWidget(0,i);
        xData[i] = sp->value();

        sp = (QDoubleSpinBox*)gridTable->cellWidget(1,i);
        zData[i] = sp->value();
    }

 nx =  nx_;
 nz = nz_;
 L = xData.last();
 d = zData.first();
 sz = sz_;

}

void customGrid::showGrid(){
    computeZ();
    QCustomPlot *cPlot = new QCustomPlot;
    QWidget *plotWindow = new QWidget;
    QHBoxLayout *plotWindow_layout = new QHBoxLayout;
    plotWindow_layout->addWidget(cPlot);
    plotWindow->setLayout(plotWindow_layout);

    plotWindow->resize(800,400);


    plotWindow->show();
    cPlot->clearGraphs();
    QVector<double> zTmp;zTmp.resize(nx);

        for (int k=0;k<nz;k++){

            for (int i=0;i<nx;i++){
                zTmp[i]=z[i*nz+k];

            }
            cPlot->addGraph();
            cPlot->graph()->setData(x, zTmp);
            cPlot->graph()->setLineStyle(QCPGraph::lsLine);
            cPlot->graph()->setScatterStyle(QCPScatterStyle::ssCross);

        }

        cPlot->xAxis->setLabel("x [m]");
        cPlot->yAxis->setLabel("z [m]");
        cPlot->xAxis->setRange(0, xData.last()+10);
        cPlot->rescaleAxes();


    cPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    cPlot->replot();

}

void customGrid::showGrid(double h, double L,int nx, int nz, bool sz){


    QCustomPlot *cPlot = new QCustomPlot;
    QWidget *plotWindow = new QWidget;
    QHBoxLayout *plotWindow_layout = new QHBoxLayout;
    plotWindow_layout->addWidget(cPlot);
    plotWindow->setLayout(plotWindow_layout);

    plotWindow->resize(800,400);
    //    plotWindow->setMinimumHeig/ht(600);


    plotWindow->show();
    cPlot->clearGraphs();


        double dx = L/(nx-1);
        std::vector<double> z;
        z.resize(nz);
        double K;
        if (sz){
            for (int k=0;k<nz;k++) {
                K =k;
                z[k] =-1*h+h*sin(K/(nz-1)*0.5*3.1415);

            }
        }else {
            for (int k=0;k<nz;k++) {
                K =k;
                z[k] =-1*h+h*K/(nz-1);

            }

        }

        QVector<double> x(nx), y(nx);


        for (int k=0;k<nz;k++){

            for (int i=0;i<nx;i++){
                x[i] = i*dx;
                y[i] = z[k];
            }
            cPlot->addGraph();
            cPlot->graph()->setData(x, y);
            cPlot->graph()->setLineStyle(QCPGraph::lsLine);
            cPlot->graph()->setScatterStyle(QCPScatterStyle::ssCross);

        }

        cPlot->xAxis->setLabel("x [m]");
        cPlot->yAxis->setLabel("z [m]");
        cPlot->xAxis->setRange(0, L+10);
        cPlot->yAxis->setRange(1, -1*h-1);


    cPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    cPlot->replot();

}

double customGrid::depth(double x){
    double h(999999999);
if (x<=xData.first()){
    h=zData.first();
} else if (x>=xData.last()) {
    h = zData.last();
} else {
    for (int i=0;i<zData.size()-1;i++){

        if ((xData[i+1]>x)&(xData[i]<x)){
            h = (zData[i+1]-zData[i])/(xData[i+1]-xData[i])*(x-xData[i])+zData[i];
        } else if (xData[i+1]==x){
            h = zData[i+1];

        }
    }
}

    return h;
}

void customGrid::generateGrid(){
    x.clear();x.resize(nx);
    z.clear();z.resize(nx*nz);
    h.clear();h.resize(nx);

    double dx = L/(nx-1);

    for (int i = 0;i<nx;i++){
        x[i] = i*dx;
        h[i] = depth(x[i]);
    }
    computeGradients(x,h);
}

void customGrid::computeZ(){
    double K;
    if (sz){
        for (int i = 0;i<nx;i++){

            for (int k =0;k<nz;k++){
                K = k;
                z[i*nz+k] = -1*h[i]+h[i]*sin(K/(nz-1)*0.5*3.1415);
            }
        }
    } else {
        for (int i = 0;i<nx;i++){

            for (int k =0;k<nz;k++){
                K = k;
                z[i*nz+k] = -1*h[i]+h[i]*K/(nz-1);
            }
        }

    }
}


void customGrid::computeGradients(QVector<double> x_,QVector<double> h_){
    hx.clear();hx.resize(nx);
    hxx.clear();hxx.resize(nx);

    for (int i=1;i<nx-1;i++){
        hx[i] = (h_[i+1] - h_[i-1])/(x_[i+1]-x_[i-1]);
        hxx[i]= (h_[i+1] -2*h_[i] + h_[i-1])/pow((x_[i+1]-x_[i]),2);

    }
    hx[0] = hx[1]; hx[nz-1]= hx[nz-1];
    hxx[0] = hxx[1]; hxx[nz-1]= hxx[nz-2];

}

void customGrid::laplaceSmoothing(){

   int N=5; // For simplicity we hardcode N, then the user don't have to make a choice
   double tmp=0;
   for (int smoothSteps=0;smoothSteps<5;smoothSteps++){
       for (int i=N/2;i<nx-N/2;i++){// we only smooth the inner part of the grid

           for (int j=-N/2;j<N/2+1;j++){
               tmp+=h[i+j];
           }
           h[i] = tmp/N;

           tmp = 0;
       }
   }
computeGradients(x,h); // After smoothing we update the gradients
}

void customGrid::writeGridToFile(QString fileName){

    QFile myfile(fileName.simplified());
    myfile.open (QIODevice::WriteOnly);
    bool notOpen = myfile.isOpen();
    bool norWrite = myfile.isWritable();
    if(!notOpen || !norWrite ){
        std::cout << "something is wrong" << std::endl;
    }
    QTextStream outStream(&myfile);
    outStream << "Custom bathymetri (2D), generated by OCW3D (GUI)\n";
    outStream << 0 << "\n";
    outStream << h.first() << " " << hx.first() << " " << hxx.first()
              << " 0 0 \n"; // ghost grid
    for (int i=0;i<nx;i++){
        outStream << h[i] << " " << hx[i] << " " << hxx[i] << " 0 0 \n";
    }
       outStream << h.last() << " " << hx.last() << " " << hxx.last() << " 0 0 \n"; // ghost grid
    myfile.close();


    // Write a file which is easy to plot
    QFile myfile2("bathymetry.plot");
    myfile2.open (QIODevice::WriteOnly);
    bool notOpen2 = myfile2.isOpen();
    bool norWrite2 = myfile2.isWritable();
    if(!notOpen2 || !norWrite2 ){
        std::cout << "something is wrong" << std::endl;
    }
    QTextStream outStream2(&myfile2);
    outStream2 << "Custom bathymetri (2D), generated by OCW3D (GUI): For plotting only \n";

    for (int i=0;i<nx;i++){
        outStream2 << x[i] << " " << h[i]  << "\n";
    }
    myfile2.close();



}
