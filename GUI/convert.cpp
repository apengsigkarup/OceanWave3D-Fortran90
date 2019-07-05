    #include "convert.h"
#include "math.h"
#include "qprogressbar.h"

// Constructor
convert::convert()
{

}
// clear memory
void convert::clearMem(){

    fileName = ' ';
    xbeg = 0;
    xend = 0;
    xstride = 0;
    ybeg = 0;
    yend = 0;
    ystride = 0;
    tbeg = 0;
    tend = 0;
    tstride = 0;
    dt = 0;
    nz = 0;
    nx = 0;
    ny = 0;
    nt = 0;

   std::vector<double> emptyVector;

//    x.swap(emptyVector);
//    y.swap(emptyVector);
//    h.swap(emptyVector);
//    hy.swap(emptyVector);
//    hx.swap(emptyVector);
//    sigma.swap(emptyVector);
//    t.swap(emptyVector);
//    F.swap(emptyVector);
//    eta.swap(emptyVector);
//    etax.swap(emptyVector);
//    etay.swap(emptyVector);
//    phi.swap(emptyVector);
//    u.swap(emptyVector);
//    v.swap(emptyVector);
//    w.swap(emptyVector);
//    uz.swap(emptyVector);
//    vz.swap(emptyVector);
   std::vector<double>().swap(x);
   std::vector<double>().swap(y);
   std::vector<double>().swap(h);
   std::vector<double>().swap(hy);
   std::vector<double>().swap(hx);
   std::vector<double>().swap(sigma);
   std::vector<double>().swap(t);
   std::vector<double>().swap(F);
   std::vector<double>().swap(eta);
   std::vector<double>().swap(etax);
   std::vector<double>().swap(etay);
   std::vector<double>().swap(phi);
   std::vector<double>().swap(u);
   std::vector<double>().swap(v);
   std::vector<double>().swap(w);
   std::vector<double>().swap(uz);
   std::vector<double>().swap(vz);

}
// Destructor
convert::~convert(){
}


void convert::read(QString file, QProgressBar* Progressbar)
{
    fileName = file;
    // Open fstream
    fileStream.open(fileName.toLatin1().data(), std::ios::binary | std::ios::in);

    if (!fileStream) // Verify that the file was open successfully
    {
        QMessageBox msgBox;
        msgBox.setText("Unknown input file");
        msgBox.exec();
        return;
    }
    fileStream.read((char *)&junk,sizeof(junk));
    fileStream.read((char *)&xbeg,sizeof(xbeg));
    fileStream.read((char *)&xend,sizeof(xend));
    fileStream.read((char *)&xstride,sizeof(xstride));
    fileStream.read((char *)&ybeg,sizeof(ybeg));
    fileStream.read((char *)&yend,sizeof(yend));
    fileStream.read((char *)&ystride,sizeof(ystride));
    fileStream.read((char *)&tbeg,sizeof(tbeg));
    fileStream.read((char *)&tend,sizeof(tend));
    fileStream.read((char *)&tstride,sizeof(tstride));
    fileStream.read((char *)&dt,sizeof(dt));
    fileStream.read((char *)&nz,sizeof(nz));

    fileStream.read((char *)&junk,2*sizeof(junk));

    nx = floor((xend-xbeg)/xstride)+1;
    ny = floor((yend-ybeg)/ystride)+1;
    nt = floor((tend-tbeg)/tstride)+1;

    //    tmp = new double[nx*ny*std::max(nz,5)];

//    x = new double[nx*ny];
    x.resize(nx*ny);
    y.resize(nx*ny);
    h.resize(nx*ny);
    hx.resize(nx*ny);
    hy.resize(nx*ny);
    sigma.resize(nz);

    for (int i=0;i<nx*ny;i++){
        fileStream.read((char *)&x[i],sizeof(dt));
        fileStream.read((char *)&y[i],sizeof(dt));
        fileStream.read((char *)&h[i],sizeof(dt));
        fileStream.read((char *)&hx[i],sizeof(dt));
        fileStream.read((char *)&hy[i],sizeof(dt));

    }


    fileStream.read((char *)&junk,2*sizeof(junk));

    // Read sigma
    for (int k=0;k<nz;k++){
        fileStream.read((char *)&sigma[k],sizeof(dt));
    }
    fileStream.read((char *)&junk,2*sizeof(junk));

    // Initialize fields
    t.resize(nt);
    n_xy  = nx*ny*sizeof(dt);
    n_xyz = nx*ny*nz*sizeof(dt);
    tmp_xy = new double[nx*ny];
    tmp_xyz = new double[nx*ny*nz];
    eta.resize(nt*nx*ny);
    etax.resize(nt*nx*ny);
    etay.resize(nt*nx*ny);
    phi.resize(nt*nx*ny*nz);
    u.resize(nt*nx*ny*nz);
    v.resize(nt*nx*ny*nz);
    w.resize(nt*nx*ny*nz);
    uz.resize(nt*nx*ny*nz);
    vz.resize(nt*nx*ny*nz);


int updateInterval = nt/200>1 ? nt/200 : 1;

    for (int i=0;i<nt;i++){

        if (i % updateInterval == 0){
            Progressbar->setValue((i+2.0)/nt*100);

        }
        // read eta
        fileStream.read((char *)tmp_xy,n_xy);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny;j++){eta[j+i*nx*ny]=tmp_xy[j];}


        //         read etax
        fileStream.read((char *)tmp_xy,n_xy);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny;j++){etax[j+i*nx*ny]=tmp_xy[j];}

        // read etay
        fileStream.read((char *)tmp_xy,n_xy);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny;j++){etay[j+i*nx*ny]=tmp_xy[j];}

        // read phi
        fileStream.read((char *)tmp_xyz,n_xyz);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny*nz;j++){phi[j+i*nx*ny*nz]=tmp_xyz[j];}

        // read u
         fileStream.read((char *)tmp_xyz,n_xyz);
         fileStream.read((char *)&junk,2*sizeof(junk));
         for (int j=0;j<nx*ny*nz;j++){u[j+i*nx*ny*nz]=tmp_xyz[j];}

        // read v
        fileStream.read((char *)tmp_xyz,n_xyz);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny*nz;j++){v[j+i*nx*ny*nz]=tmp_xyz[j];}

        // read w
        fileStream.read((char *)tmp_xyz,n_xyz);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny*nz;j++){w[j+i*nx*ny*nz]=tmp_xyz[j];}

        // read uz
        fileStream.read((char *)tmp_xyz,n_xyz);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny*nz;j++){uz[j+i*nx*ny*nz]=tmp_xyz[j];}

        // read vz
        fileStream.read((char *)tmp_xyz,n_xyz);
        fileStream.read((char *)&junk,2*sizeof(junk));
        for (int j=0;j<nx*ny*nz;j++){vz[j+i*nx*ny*nz]=tmp_xyz[j];}

        t[i] = (i+1)*dt;
    }

    fileStream.close();


}

void convert::netCfd(){



}


int convert::matlab(){
#if MATLAB>0
    MATFile *pmat;
    mxArray *time_mat, *eta_mat, *x_mat,*y_mat,*u_mat,*v_mat,*w_mat,*uz_mat,*vz_mat,*sigma_mat,*h_mat;

    const char *file = "OCW3D.mat";
    //char str[BUFSIZE];
    int status;




    // Open MAT file
    printf("Creating file %s...\n\n", file);
    pmat = matOpen(file, "w");
    if (pmat == NULL) {
        printf("Error creating file %s\n", file);
        printf("(Do you have write permission in this directory?)\n");
        return(EXIT_FAILURE);
    }

    // TIME
    time_mat = mxCreateDoubleMatrix(nt,1,mxREAL);
    if (time_mat == NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(time_mat)),t,nt*sizeof(t));
    status = matPutVariable(pmat, "t", time_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    //// eta
    eta_mat = mxCreateDoubleMatrix(nx*ny,nt,mxREAL);
    if (eta_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(eta_mat)), (void *)eta,nt*nx*ny*sizeof(eta));
    status = matPutVariable(pmat, "eta", eta_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    // h
    h_mat = mxCreateDoubleMatrix(nx*ny,1,mxREAL);
    if (h_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(h_mat)), (void *)h,nx*ny*sizeof(h));
    status = matPutVariable(pmat, "h", h_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    // sigma
    sigma_mat = mxCreateDoubleMatrix(nz,1,mxREAL);
    if (sigma_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(sigma_mat)), (void *)sigma,nz*sizeof(sigma));
    status = matPutVariable(pmat, "sigma", sigma_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    // x
    x_mat = mxCreateDoubleMatrix(nx,1,mxREAL);
    if (x_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(x_mat)), (void *)x,nx*sizeof(x));
    status = matPutVariable(pmat, "x", x_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

    // y
    y_mat = mxCreateDoubleMatrix(ny,1,mxREAL);
    if (y_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    memcpy((void *)(mxGetPr(y_mat)), (void *)y,ny*sizeof(y));
    status = matPutVariable(pmat, "y", y_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }
    // u
    u_mat = mxCreateDoubleMatrix(nx*ny*nz,nt,mxREAL);
    if (u_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    double *pointer; pointer = mxGetPr(u_mat);
    for (int index=0;index<nt*nx*ny*nz;index++){ pointer[index] = u[index];}
    status = matPutVariable(pmat, "u", u_mat);

    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

// v
    v_mat = mxCreateDoubleMatrix(nx*ny*nz,nt,mxREAL);
    if (v_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    pointer = mxGetPr(v_mat); for (int index=0;index<nt*nx*ny*nz;index++){ pointer[index] = v[index];}
    status = matPutVariable(pmat, "v", v_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

// w
    w_mat = mxCreateDoubleMatrix(nx*ny*nz,nt,mxREAL);
    if (w_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    pointer = mxGetPr(w_mat); for (int index=0;index<nt*nx*ny*nz;index++){ pointer[index] = w[index];}
    status = matPutVariable(pmat, "w", w_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

// uz
    uz_mat = mxCreateDoubleMatrix(nx*ny*nz,nt,mxREAL);
    if (uz_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    pointer = mxGetPr(uz_mat); for (int index=0;index<nt*nx*ny*nz;index++){pointer[index] = uz[index];
    }
    status = matPutVariable(pmat, "uz", uz_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }

// vz
    vz_mat = mxCreateDoubleMatrix(nx*ny*nz,nt,mxREAL);
    if (vz_mat== NULL) {
        printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
        printf("Unable to create mxArray.\n");
        return(EXIT_FAILURE);
    }
    pointer = mxGetPr(vz_mat); for (int index=0;index<nt*nx*ny*nz;index++){ pointer[index] = vz[index];
    }
    status = matPutVariable(pmat, "vz", vz_mat);
    if (status != 0) {
        printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
        return(EXIT_FAILURE);
    }


    /* clean up */
    mxDestroyArray(time_mat);
    mxDestroyArray(eta_mat);
    mxDestroyArray(x_mat);
    mxDestroyArray(y_mat);
    mxDestroyArray(u_mat);
    mxDestroyArray(v_mat);
    mxDestroyArray(w_mat);
    mxDestroyArray(uz_mat);
    mxDestroyArray(vz_mat);
    mxDestroyArray(sigma_mat);
    mxDestroyArray(h_mat);


    if (matClose(pmat) != 0) {
        printf("Error closing file %s\n",file);
        return(EXIT_FAILURE);
    }
    return 1;
#else
    return 0;
#endif

}
void convert::gradient(std::vector<double> argout,const std::vector<double> argin,const double dt,const int n) const
{

    // First and last index (one-sided)
    //
    argout[0] = (argin[1]-argin[0])/dt;
    argout[n-1] = (argin[n-1]-argin[n-2])/dt;

    for (int j=1;j<n-2;j++){
        argout[j] = (argin[j+1]-argin[j-1])/(2*dt);
    }


}

void convert::gradient(double argout[],const std::vector<double> argin,const double dt,const int n) const
{

    // First and last index (one-sided)
    //
    argout[0] = (argin[1]-argin[0])/dt;
    argout[n-1] = (argin[n-1]-argin[n-2])/dt;

    for (int j=1;j<n-2;j++){
        argout[j] = (argin[j+1]-argin[j-1])/(2*dt);
    }


}
void convert::gradient(double argout[],const double argin[],const double dt,const int n) const
{

    // First and last index (one-sided)
    //
    argout[0] = (argin[1]-argin[0])/dt;
    argout[n-1] = (argin[n-1]-argin[n-2])/dt;

    for (int j=1;j<n-2;j++){
        argout[j] = (argin[j+1]-argin[j-1])/(2*dt);
    }


}

void convert::force(int nn,double D,double rho,double Cd,double Cm){
    // Define Pi
    //
    const double PI = 4.0*atan(1.0);


    // To cover-up for previous lazyness are we now making local
    // copies of the fields we need - sorry!
    //
    std::vector<double> tmp_;
    tmp_.resize(nt);
    for (int i=0;i<nt;i++){
        tmp_[i] = eta[i*nx*ny+nn];
    }

    // deta/dt
    //
    std::vector<double> detadt;
   detadt.resize(nt);
   gradient(detadt,tmp_,dt,nt);

    // dz/dt
    //
   Double2d dzdt(boost::extents[nz][nt]);
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
         dzdt[k][i] = detadt[i]*sigma[k];
        }
    }

    // du/dt on the sigma grid
    //
    Double2d dudt(boost::extents[nz][nt]);
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
            tmp_[i] = u[i*nx*ny*nz+nn*nz+k];
        }
        gradient(&dudt[k][0],tmp_,dt,nt); // To-do: this is super sluppy and error prone
     }



    // Acceleration
    //
    Double2d acc(boost::extents[nz][nt]);
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
            acc[k][i] = dudt[k][i]-uz[i*nx*ny*nz+nn*nz+k]*dzdt[k][i];
        }
    }

    // The inline force
    //
    double dz,Fd,Fi;
    int index;
    F.resize(nt);

    for (int i=0;i<nt;i++){
        F[i] = 0;
        for (int k=1;k<nz-1;k++){// we ignore the ghost points
            index = i*nx*ny*nz+nn*nz+k;
            dz = (sigma[k+1]-sigma[k])*(eta[i*nx*ny+nn]+h[nn]);

            Fd = 0.5*rho*Cd*D*(u[index]*std::abs(u[index])+u[index+1]*std::abs(u[index+1]))/2;
            Fi = rho*Cm*(PI/4)*pow(D,2)*(acc[k][i]+acc[k+1][i])/2;

            F[i] += (Fd+Fi)*dz;
        }

    }


    QFileInfo fileInfo =  QFileInfo(fileName);
    std::ofstream morisonForce;
    morisonForce.open(fileInfo.baseName().toLatin1() + ".force");
    morisonForce << "time F" << std::endl;
    QVector<double> QV_t,QV_F;
    QV_t.resize(nt);
    QV_F.resize(nt);
    for (int i=1;i<nt-1;i++){
        morisonForce << std::setiosflags(std::ios::fixed) << std::setprecision(10) << t[i] << " " << F[i] << std::endl;
        QV_t[i] = t[i];
        QV_F[i] = F[i];
    }
    morisonForce.close();

    // plot solution
    QCustomPlot *cPlot = new QCustomPlot;
    QWidget *plotWindow = new QWidget;
    QHBoxLayout *plotWindow_layout = new QHBoxLayout;
    plotWindow_layout->addWidget(cPlot);
    plotWindow->setLayout(plotWindow_layout);

    plotWindow->resize(800,400);

    plotWindow->show();
    cPlot->clearGraphs();

    cPlot->addGraph();
    cPlot->graph()->setData(QV_t,QV_F);
    cPlot->xAxis->setLabel("Time, t s");
    cPlot->yAxis->setLabel("Inline force, F N");
    QVector<double>::iterator Xmin = std::min_element(QV_t.begin(), QV_t.end());
    QVector<double>::iterator Xmax = std::max_element(QV_t.begin(), QV_t.end());
    QVector<double>::iterator Ymin = std::min_element(QV_F.begin(), QV_F.end());
    QVector<double>::iterator Ymax = std::max_element(QV_F.begin(), QV_F.end());
    cPlot->xAxis->setRange(*Xmin,*Xmax);
    cPlot->yAxis->setRange(*Ymin,*Ymax);
    cPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    cPlot->replot();


}


void convert::ascii(QString filepath,int location){
QFileInfo fileInfo =  QFileInfo(filepath);


std::ofstream oStream;
    oStream.open(fileInfo.baseName().toLatin1() + ".eta");
    if (oStream.is_open()){
        if (location < 0){
            // write the header lines
            oStream << "t x[t,1] x[t,2] x[t,3]: Second row is x-locations: Third row is water depths\n";
//            oStream << "x = ";
            for (int j=0;j<nx*ny;j++){oStream << x[j] << " ";}oStream << "\n";
//            oStream << "h = ";
            for (int j=0;j<nx*ny;j++){oStream << h[j] << " ";}oStream << "\n";


            for (int i=0;i<nt-1;i++){
                oStream << t[i] << " ";

                for (int j=0;j<nx*ny;j++){

                    oStream << eta[j+i*nx*ny] << " ";

                }
                oStream << "\n";
            }
        } else{

            // write the header lines
            oStream << "t x[t] \n ";
            oStream << "x = " << x[location] << "\n";
            oStream << "h = " << h[location] << "\n";


            for (int i=0;i<nt-1;i++){
                oStream << t[i] << " " << eta[location+i*nx*ny] << "\n";
            }


        }
        oStream.close();
    }
}


