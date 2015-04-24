#include "convert.h"
#include "qstring.h"
#include "math.h"
#include "qprogressbar.h"
//#include "mat.h"

// Constructor
convert::convert()
{



}

// Destructor
convert::~convert(){
}


void convert::read(QString file, QProgressBar* Progressbar)
{
    fileName = file;
    Progressbar->setVisible(true);
    Progressbar->setValue(0);
    // Open fstream
    fileStream.open(fileName.toLatin1().data(), std::ios::binary | std::ios::in);

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

    x = new double[nx*ny];
    y = new double[nx*ny];
    h = new double[nx*ny];
    hx = new double[nx*ny];
    hy = new double[nx*ny];
    sigma = new double[nz];

    for (int i=0;i<nx*ny;i++){
        fileStream.read((char *)&x[i],sizeof(dt));
        fileStream.read((char *)&y[i],sizeof(dt));
        fileStream.read((char *)&h[i],sizeof(dt));
        fileStream.read((char *)&hx[i],sizeof(dt));
        fileStream.read((char *)&hy[i],sizeof(dt));

    }


    fileStream.read((char *)&junk,2*sizeof(junk));

    // Read sigma
    fileStream.read((char *)sigma,nz*sizeof(dt));
    fileStream.read((char *)&junk,2*sizeof(junk));

    // Initialize fields
    t = new double[nt];
    n_xy  = nx*ny*sizeof(dt);
    n_xyz = nx*ny*nz*sizeof(dt);
    tmp_xy = new double[nx*ny];
    tmp_xyz = new double[nx*ny*nz];
    eta   = new double[nt*nx*ny];
    etax  = new double[nt*nx*ny];
    etay  = new double[nt*nx*ny];
    phi   = new double[nt*nx*ny*nz];
    u     = new double[nt*nx*ny*nz];
    v     = new double[nt*nx*ny*nz];
    w     = new double[nt*nx*ny*nz];
    uz    = new double[nt*nx*ny*nz];
    vz    = new double[nt*nx*ny*nz];



    for (int i=0;i<nt;i++){

        Progressbar->setValue((i+1)/nt*100);

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

        t[i] = i*dt;
    }

    fileStream.close();


}

void convert::netCfd(){



}


int convert::matlab(){
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


}

void convert::gradient(double argout[],const double argin[],const double dt,const int n){

    // First and last index (one-sided)
    //
    argout[0] = (argin[1]-argin[0])/dt;
    argout[n-1] = (argin[n-1]-argin[n-2])/dt;

    for (int j=1;j<n-2;j++){
        argout[j] = (argin[j+1]-argin[j-1])/(2*dt);
    }


}
void convert::force(double x0,double D,double rho,double Cd,double Cm){
    // Define Pi
    //
    const double PI = 4.0*atan(1.0);
    double dx = x[1]-x[0];
    double stencil[3];
    stencil[0] = 1/(dx*dx);stencil[1]=-2/(dx*dx);stencil[2]=stencil[0];
    double A = PI*pow((D/2),2);

    // Find nearest neighbour to x0
    //
    int nn = round((x0-xbeg*dx)/dx)+1;

    // To cover-up for previous lazyness are we now making local
    // copies of he fields we need - sorry!
    //
    double tmp_[nt];
    for (int i=0;i<nt;i++){
        tmp_[i] = eta[i*nx*ny+nn];

    }

    // deta/dt
    //
    double detadt[nt];
    gradient(detadt,tmp_,dt,nt);

    // dz/dt
    //
    double dzdt[nz][nt];
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
            dzdt[k][i] = detadt[i]*sigma[k];
        }
    }

    // du/dt on the sigma grid
    //
    double dudt[nz][nt];
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
            tmp_[i] = u[i*nx*ny*nz+nn*nz+k];
        }
        gradient(&dudt[k][1],tmp_,dt,nt); // To-do: this is super sluppy and error prone
     }



    // Acceleration
    //
    double acc[nz][nt];
    for (int k=0;k<nz;k++){
        for (int i=0;i<nt;i++){
            acc[k][i] = dudt[k][i]-uz[i*nx*ny*nz+nn*nz+k]*dzdt[k][i];
        }
    }



    // The inline force
    //
    double dz,Fd,Fi;
    int index;
    F = new double[nt];

    for (int i=0;i<nt;i++){
        F[i] = 0;
        for (int k=1;k<nz-1;k++){// we ignore thte ghost point
            index = i*nx*ny*nz+nn*nz+k;
            dz = (sigma[k+1]-sigma[k])*(eta[i*nx*ny+nn]+h[nn]);

            Fd = 0.5*rho*Cd*D*(u[index]*std::abs(u[index])+u[index+1]*std::abs(u[index+1]))/2;
            Fi = rho*Cm*(PI/4)*pow(D,2)*(acc[k][i]+acc[k+1][i])/2;

            F[i] += (Fd+Fi)*dz;
            std::cout << Fd << " " << index <<  std::endl;
        }

    }
    std::cout << "x0 = " << x[nn] << " " << nt*nx*nz<< std::endl;
    std::ofstream debug;
    debug.open("debug.force");
    for (int i=0;i<nt;i++){
        debug << std::setiosflags(std::ios::fixed) << std::setprecision(10) << t[i] << " " << F[i] << std::endl;
    }
    debug.close();



}


void convert::ascii(){
    std::cout << "convertign to ascii" << std::endl;
    std::ofstream oStream;

    oStream.open("Kinematics01.eta");
    if (oStream.is_open()){

        // write the header lines
        oStream << "t x[t,1] x[t,2] x[t,3]\n";
        oStream << " ";
        for (int j=0;j<nx*ny;j++){oStream << x[j] << " ";}oStream << "\n";


        for (int i=0;i<nt;i++){
            oStream << t[i] << " ";

            for (int j=0;j<nt;j++){

                oStream << eta[j+i*nx*ny] << " ";

            }
            oStream << "\n";
        }
        oStream.close();
    }else{

        std::cout << "Error open file" << std::endl;
    }
}
