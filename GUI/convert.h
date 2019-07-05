#ifndef CONVERT_H
#define CONVERT_H

#include "qstring.h"
#include "versions.h"
#include <QString>
#include <fstream>
#include <iterator>
#include <iostream>
#include <QProgressBar>
#include <iomanip>
#include <QFileInfo>
#include "math.h"
#include <cmath>
#include "boost/multi_array.hpp"
#include <vector>
#include "QTwidgets/qcustomplot.h"
#include <QHBoxLayout>
#include <QMessageBox>

#if MATLAB>0
    #include "mat.h"
    #include "matrix.h"
#endif

typedef boost::multi_array<double, 2> Double2d;
typedef boost::multi_array<double, 3> Double3d;

class convert
{
public:
    convert();
    ~convert();
    void read(QString,QProgressBar *);
    void netCfd();
    void force(int,double,double,double,double);
    int matlab();
    void clearMem();
    void ascii(QString, int location);
    void gradient(double*,const double*,const double,const int) const;
    void gradient(double*,const std::vector<double>,const double,const int) const;
    void gradient(std::vector<double>,const std::vector<double>,const double,const int) const;
    QString fileName;
    int xbeg;int xend;int xstride;
    int ybeg;int yend;int ystride;
    int tbeg;int tend;int tstride;double dt;
    int nz;int nx;int ny; int nt;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> h;
    std::vector<double> hy;
    std::vector<double> hx;
    std::vector<double> sigma;
    std::vector<double> t;
    std::vector<double> F;
    std::vector<double> eta;
    std::vector<double> etax;
    std::vector<double> etay;
    std::vector<double> phi;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> uz;
    std::vector<double> vz;


private:

    // private member functions



    std::size_t n_xy;
    std::size_t n_xyz;
    std::size_t n_z;
    double* tmp_z;
    double* tmp_xy;
    double* tmp_xyz;

   int junk;
//   double* tmp;
   std::streampos pos;
   std::ifstream fileStream;
};

#endif // CONVERT_H
