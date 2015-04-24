#ifndef CONVERT_H
#define CONVERT_H

#include "qstring.h"
#include <QString>
#include <fstream>
#include <iterator>
#include <iostream>
#include <QProgressBar>
#include <iomanip>

#include "mat.h"
#include "math.h"
#include <cmath>
#include "matrix.h"

class convert
{
public:
    convert();
    ~convert();
    void read(QString,QProgressBar *);
    void netCfd();
    void force(double,double,double,double,double);
    int matlab();
    void ascii();

    QString fileName;
    int xbeg;int xend;int xstride;
    int ybeg;int yend;int ystride;
    int tbeg;int tend;int tstride;double dt;
    int nz;int nx;int ny; int nt;

    double *x;double* y;double* h;double* hy;double* hx;
    double* sigma; double* t; double* F;
    double* eta;double* etax;double* etay;double* phi;double *u;
    double* v; double* w; double* uz; double* vz;


private:

    // private member functions
    void gradient(double*,const double*,const double,const int);


    std::size_t n_xy;
    std::size_t n_xyz;
    double* tmp_xy;
    double* tmp_xyz;

   int junk;
//   double* tmp;
   std::streampos pos;
   std::ifstream fileStream;
};

#endif // CONVERT_H
