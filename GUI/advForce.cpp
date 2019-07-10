#include "advForce.h"

advForce::~advForce()
{
    delete ui;
}

advForce::advForce(const convert &data_, QWidget *parent) :
    QWidget(parent)
{
    ui = new QWidget;
    mainLayout = new QGridLayout;
    data = &data_;

    createLeftSide();
    createRightSide();

    mainLayout->addWidget(leftView2,0,0,55,100);
    mainLayout->addWidget(leftView,55,0,10,100);
    mainLayout->addWidget(leftView3,67,0,2,100);

    mainLayout->addWidget(rightView,0,250,80,250);
    ui->setLayout(mainLayout);
    ui->resize(1500,1000);
    ui->show();

    connect(nRows,SIGNAL(valueChanged(int)),this,SLOT(updateTableSize_()));
    connect(selectadvForceFile,SIGNAL(clicked()),this,SLOT(setadvForceFile_()));
    connect(selectGeometryFile,SIGNAL(clicked()),this,SLOT(setGeometryFile_()));
    connect(updatePlotAndFile,SIGNAL(clicked()),this,SLOT(plotMonoPile_()));
    connect(save,SIGNAL(clicked()),this,SLOT(save_()));
    connect(plot,SIGNAL(clicked()),this,SLOT(plot_()));
    connect(gridType,SIGNAL(currentIndexChanged(int)),this,SLOT(gridTypeChanged_()));
}

void advForce::setadvForceFile_(){


    advForceFile_fileInfo = selectFile();
    advForceFile->setText(advForceFile_fileInfo.fileName());

}

void advForce::setGeometryFile_()
{

    geometryFile_fileInfo = selectFile();
    geometryFile->setText(geometryFile_fileInfo.fileName());
    readGeometryFile();
}

QFileInfo advForce::selectFile()
{
    QString fileName = QFileDialog::getOpenFileName(
                this,
                "Select file"
                );

    QFile inputFile(fileName);

    return inputFile;

}

void advForce::createRightSide()
{
    rightView = new QGroupBox(tr("Geometry"));
    monopilePlot = new QCustomPlot(this);
    QGridLayout *layout = new QGridLayout;

    monopilePlot->addGraph(monopilePlot->xAxis,monopilePlot->yAxis);


    layout->addWidget(monopilePlot,0,0,80,250);

    rightView->setLayout(layout);


}

void advForce::createLeftSide()
{
    // initialize
    leftView                = new QGroupBox(tr("Output settings"));
    leftView2               = new QGroupBox(tr("Geometry settings"));
    leftView3               = new QGroupBox(tr("Save and plot"));

    QGridLayout *layout     = new QGridLayout;
    QGridLayout *layout2    = new QGridLayout;
    QHBoxLayout *layout3    = new QHBoxLayout;

    advForceFile               = new QLineEdit(this);
    selectadvForceFile         = new QPushButton(this);
    selectadvForceFile_lable   = new QLabel(this);
    nradvForceLevels_label     = new QLabel(tr("Number of levels in constant grid"));
    xLocation               = new QComboBox(this);
    xLocation_label         = new QLabel(tr("x-location"));
    nradvForceLevels           = new QSpinBox(this);
    gridType                = new QComboBox(this);



    selectGeometryFile_label= new QLabel(tr("GeometryFile (*.monopile)"));
    geometryFile            = new QLineEdit(this);
    geometryHeader          = new QLineEdit(this);
    geometryHeader_label    = new QLabel("Description");
    selectGeometryFile      = new QPushButton(tr("select"));
    forceCoeffsAndDiameter  = new QTableWidget(this);
    nRows                   = new QSpinBox(this);
    nVerticalLayers_label   = new QLabel(tr("Number of vertical levels (geometry)"));
    density_label           = new QLabel(tr("Density kg/m^3"));
    density                 = new QDoubleSpinBox(this);
    updatePlotAndFile       = new QPushButton(tr("plot and save geometry"));
    save                    = new QPushButton(tr("Save"));
    plot                    = new QPushButton(tr("Plot"));


    // settings
    layout->setColumnMinimumWidth(0,1);
    layout->setColumnStretch(1,10);
    layout2->setColumnMinimumWidth(0,1);

    selectadvForceFile->setText(tr("select"));
    selectadvForceFile_lable->setText("Output file");
    nradvForceLevels->setValue(50);
    nradvForceLevels_label->setVisible(false); // change to falase
    nradvForceLevels->setVisible(false); // change to false
    forceCoeffsAndDiameter->setGeometry(QRect(0,0,15,20));
    nRows->setValue(1);
    updateTableSize(forceCoeffsAndDiameter,nRows->value(),4);
    forceCoeffsAndDiameter->setHorizontalHeaderLabels(QString("z-level;Diameter;Cd;Cm").split(";"));
    density->setMaximum(10000);
    density->setValue(1025);
    updatePlotAndFile->setFont(QFont("Arial", 12, QFont::Bold));
    save->setFont(QFont("Arial", 12, QFont::Bold));
    plot->setFont(QFont("Arial", 12, QFont::Bold));

    // populate ComboBox

    if (!data->x.empty()){

        for (int i=0;i<data->nx;i++)
        {
            xLocation->addItem(QString::number(data->x[i],'f',2));
        }
    }

    gridType->addItem("Sigma grid");
    gridType->addItem("Cartesian grid (not implemented yet)");
    gridType->setFixedSize(200,24);
    // set layout

    layout->addWidget(selectadvForceFile_lable,0,0,1,1);
    layout->addWidget(advForceFile,1,0,1,27);
    layout->addWidget(selectadvForceFile,1,27,1,1);
    layout->addWidget(xLocation_label,2,0,1,1);
    layout->addWidget(xLocation,3,0,1,1);
    layout->addWidget(gridType,3,1,1,1);
    layout->addWidget(nradvForceLevels_label,2,20,1,20);
    layout->addWidget(nradvForceLevels,3,20,1,1);

    leftView->setLayout(layout);

    layout2->addWidget(selectGeometryFile_label,0,0,1,27);
    layout2->addWidget(geometryFile,1,0,1,27);
    layout2->addWidget(selectGeometryFile,1,27,1,3);
    layout2->addWidget(geometryHeader_label,4,0,1,2);
    layout2->addWidget(geometryHeader,5,0,1,27);
    layout2->addWidget(nVerticalLayers_label,6,0,1,2);
    layout2->addWidget(nRows,7,0,1,1);
    layout2->spacerItem();
    layout2->addWidget(forceCoeffsAndDiameter,6,2,22,28);
    layout2->addWidget(density_label,8,0,1,1);
    layout2->addWidget(density,9,0,1,1);
    layout2->addWidget(updatePlotAndFile,100,0,2,30);


    leftView2->setLayout(layout2);


    layout3->addWidget(save);
    layout3->addWidget(plot);

    leftView3->setLayout(layout3);

    Fout.resize(data->nt);


}

void advForce::gridTypeChanged_(){
    if (gridType->currentIndex()==1) {
        QMessageBox msgBox;
        msgBox.setText("Sorry, This is not implemented yet");
        msgBox.exec();

gridType->setCurrentIndex(0);
//        nradvForceLevels_label->setVisible(true);
//        nradvForceLevels->setVisible(true);
    } else {

        nradvForceLevels_label->setVisible(false);
        nradvForceLevels->setVisible(false);
    }


}

void advForce::plotMonoPile_()
{

    int zCol = 0;
    int diameterCol = 1;

    QVector<double> x,x2,y,y2;
    QDoubleSpinBox *sp;
    for (int i =0;i<forceCoeffsAndDiameter->rowCount();i++){
        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,diameterCol);
        x.push_back(sp->value()/2);
        x2.push_back(sp->value()/-2);

        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,zCol);
        y.push_back(sp->value());
        y2.push_back(sp->value());
    }

    // prepare for new plot
    monopilePlot->clearGraphs();

    // compute graph limits
    double ylim_low = (*std::min_element(y.begin(),y.end()));
    double ylim_up  = (*std::max_element(y.begin(),y.end()))+1;
    double xlim     = (ylim_up-ylim_low);

    QVector<double> tmpy(2,ylim_low),tmpx(2);
    tmpx[0] = -xlim;tmpx[1] = xlim;

    // add the seabed for ylim_low
    QPen pen_seaBed;
    pen_seaBed.setColor(QColor(0,0,0));
    pen_seaBed.setStyle(Qt::SolidLine);
    pen_seaBed.setWidthF(12);
    monopilePlot->addGraph();
    monopilePlot->graph()->addData(tmpx,tmpy);
    monopilePlot->graph()->setPen(pen_seaBed);

    // add still water level as dashed line
    QPen pen_swl;
    tmpy.fill(0);
    pen_swl.setColor(QColor(0,0,20));
    pen_swl.setStyle(Qt::DashLine);
    pen_swl.setWidth(1);
    monopilePlot->addGraph();
    monopilePlot->graph()->addData(tmpx,tmpy);
    monopilePlot->graph()->setPen(pen_swl);
    // Plot the monopile
    // This is not very pretty, but it is made this way due to an unexpected bhaviour of qcustomplot
    for (int i=0;i<x.size()-1;i++){
        tmpx[0] = x.at(i) ; tmpx[1] = x.at(i+1);
        tmpy[0] = y.at(i); tmpy[1] = y.at(i+1);
        monopilePlot->addGraph(monopilePlot->yAxis,monopilePlot->xAxis); // we swap the axis for the brush to work
        monopilePlot->graph()->setBrush(QBrush(QColor(0, 0, 255, 50)));
        monopilePlot->graph()->setData(tmpy,tmpx);


        tmpx[0] = x2.at(i); tmpx[1] = x2.at(i+1);
        monopilePlot->addGraph(monopilePlot->yAxis,monopilePlot->xAxis); // we swap the axis for the brush to work
        monopilePlot->graph()->setData(tmpy,tmpx);
        monopilePlot->graph()->setBrush(QBrush(QColor(0, 0, 255, 50)));
    }


    monopilePlot->xAxis->setRange(-xlim, xlim);
    monopilePlot->yAxis->setRange(ylim_low,ylim_up);
    monopilePlot->xAxis->setLabel("Horizontal, x m");
    monopilePlot->yAxis->setLabel("Vertical, z m");
    monopilePlot->replot();

    saveGeometry();
}

void advForce::updateTableSize_(){

    updateTableSize(forceCoeffsAndDiameter,nRows->value(),4);

}

void advForce::updateTableSize(QTableWidget *table,int nRow,int nCol)
{
    QDoubleSpinBox *sp;
    // check if there are too many rows
    if (table->rowCount()>nRow){
        for (int i = table->rowCount();i>nRow;i--){
            table->removeRow(i-1);
        }
    }

    // check if there are too many columns
    if (table->columnCount()>nCol){
        for (int j = table->columnCount();j>nCol;j--){
            table->removeColumn(j-1);
        }
    }

    // check if we need more rows
    if (table->rowCount()<nRow){

        for (int i = table->rowCount();i<nRow;i++){
            table->insertRow(i);

            for (int jj = 0;jj<table->columnCount();jj++){
                table->setCellWidget(i,jj,new QDoubleSpinBox(table));
                sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,jj);
                sp->setMinimum(-1000);
                sp->setMaximum(1000);
            }
        }
    }

    // check if we need more columns
    if (table->columnCount()<nCol){

        for (int j = table->columnCount();j<nCol;j++){
            table->insertColumn(j);

            for (int ii = 0;ii<table->rowCount();ii++){
                table->setCellWidget(ii,j,new QDoubleSpinBox(table));
                sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(ii,j);
                sp->setMinimum(-1000);
                sp->setMaximum(1000);
            }
        }
    }

}

void advForce::readGeometryFile()
{
    QString tmp_line;
    QStringList tmp_list;
    QDoubleSpinBox *sp;

    QFile in(geometryFile_fileInfo.filePath());
    in.open(QIODevice::ReadOnly);

    // first line is a header so we skip it (can be added as title to plot??)
    tmp_line = in.readLine();
    geometryHeader->setText(tmp_line);

    // second line is the number of lines
    tmp_line = in.readLine();
    tmp_list = tmp_line.split(" ");

    nRows->setValue(tmp_list[0].toInt());
    updateTableSize_();
    int i=0;
    while(!in.atEnd()){


        tmp_line = in.readLine();
        tmp_list = tmp_line.split(" ");

        for (int j =0;j<forceCoeffsAndDiameter->columnCount();j++){
            sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,j);
            sp->setValue(tmp_list[j].toDouble());
        }

        i++;
    }

    plotMonoPile_();
}


void advForce::computeForce(int nn){
    // Define Pi
    //
    const double PI = 4.0*atan(1.0);

    // To cover-up for previous lazyness are we now making local
    // copies of the fields we need - sorry!
    //
    std::vector<double> tmp_;
    tmp_.resize(data->nt);
    for (int i=0;i<data->nt;i++){
        tmp_[i] = data->eta[i*data->nx*data->ny+nn];
    }

    // deta/dt
    //
    std::vector<double> detadt;
    detadt.resize(data->nt);
    data->gradient(detadt,tmp_,data->dt,data->nt);

    // dz/dt
    //
    Double2d dzdt(boost::extents[data->nz][data->nt]);
    for (int k=0;k<data->nz;k++){
        for (int i=0;i<data->nt;i++){
            dzdt[k][i] = detadt[i]*data->sigma[k];
        }
    }

    // du/dt on the sigma grid
    //
    Double2d dudt(boost::extents[data->nz][data->nt]);
    for (int k=0;k<data->nz;k++){
        for (int i=0;i<data->nt;i++){
            tmp_[i] = data->u[i*data->nx*data->ny*data->nz+nn*data->nz+k];
        }
        data->gradient(&dudt[k][0],tmp_,data->dt,data->nt); // To-do: this is super sluppy and error prone
    }



    // Acceleration
    //
    Double2d acc(boost::extents[data->nz][data->nt]);
    for (int k=0;k<data->nz;k++){
        for (int i=0;i<data->nt;i++){
            acc[k][i] = dudt[k][i]-data->uz[i*data->nx*data->ny*data->nz+nn*data->nz+k]*dzdt[k][i];
        }
    }



    // The inline force
    //
    double Fd,Fi;
    std::vector<double> z_ocw;
    std::vector<double> dz;
    int index;

    F.resize(data->nz-1);
    z_ocw.resize(data->nz-1);


    // Read the table with D,cm and cd
    std::vector<double> D,Cd,Cm;
    double rho = density->value();
    readTable();

    std::ofstream morisonForce;
    morisonForce.open(advForceFile_fileInfo.fileName().toLatin1());
    morisonForce << "time [s]   Force [N]  \eta [m] \n";


    double etaMax = *std::max_element(data->eta.begin(),data->eta.end());
    std::vector<double> b_grid;
    std::vector<double> F_advForce;
    std::vector<double> dz_b;
//    if (gridType->currentIndex()==1) {
//        // we need to create a grid which is stationary in time
//        //


//        int n_advForce = nradvForceLevels->value();


//        b_grid.resize(n_advForce);
//        F_advForce.resize(n_advForce);
//        dz_b.resize(n_advForce);

//        for(int i=0;i<n_advForce;i++){
//            b_grid[i] =  std::sin(PI/2/n_advForce*i)*(data->h[nn]+etaMax+1) - data->h[nn];
//        }



//        dz_b[0] = (b_grid[1]-b_grid[0])/2; // for k=0, we integrate from z=0 to z=dz/2

//        // for all other points we
//        for (int k = 1;(unsigned)k<b_grid.size()-1;k++){
//            dz_b[k] = (b_grid[k+1]-b_grid[k])/2 + (b_grid[k]-b_grid[k-1])/2;

//        }
//        dz_b.back() =  (b_grid.back()-b_grid[b_grid.size()-2])/2;

//        morisonForce << b_grid.size() << "\n";
//        for (int k = 0;(unsigned)k<b_grid.size();k++){
//            morisonForce << std::setiosflags(std::ios::fixed) << std::setprecision(3) << b_grid[k]+data->h[nn] << " ";
//        }
//        morisonForce << "\n";
//        morisonForce << data->nt << "\n";
//    }

    for (int i=0;i<data->nt-1;i++){
        F.empty();
        z_ocw.empty();

        // For each time step we need the physical OceanWave3D grid
        //
        for (int k=1;(unsigned)k<data->sigma.size();k++){
            z_ocw[k-1] = data->sigma[k]*(data->eta[i*data->nx*data->ny+nn]+data->h[nn])-data->h[nn];
        }

        D.clear();D.resize(z_ocw.size());
        Cd.clear();Cd.resize(z_ocw.size());
        Cm.clear();Cm.resize(z_ocw.size());

        int dummy;
        interp_linear(D_table,z_table,z_ocw,D,dummy);
        interp_linear(Cd_table,z_table,z_ocw,Cd,dummy);
        interp_linear(Cm_table,z_table,z_ocw,Cm,dummy);
        // First we loop in the vertical and compute the forces on the OceanWave3D grid
        //
        for (int k=1;k<data->nz;k++){// we ignore the ghost points
            index = i*data->nx*data->ny*data->nz+nn*data->nz+k;

            Fd = 0.5*rho*Cd[k]*D[k]*data->u[index]*std::fabs(data->u[index]);
            Fi = rho*Cm[k]*(PI/4)*pow(D[k],2)*(acc[k][i]);

            F[k-1] = (Fd+Fi);//*dz[k-1];
        }



        if (gridType->currentIndex()==1) {
            // we interpolate to the advForce grid
            int id=0; // we need this later
            interp_linear(F,z_ocw,b_grid,F_advForce,id);

            // Finially we have to correct the end point
            //
            double lambda = z_ocw.back();
            double L0 = b_grid[id];
            double w = (lambda-(L0+dz_b[id]/2))/dz_b[id];
            F_advForce[id] *= (1+w);
            F_advForce[id+1] = 0;

            morisonForce << data->t[i] << " ";
            for (int k = 0;(unsigned)k<F_advForce.size();k++){
                morisonForce << std::setiosflags(std::ios::fixed) << std::setprecision(2) << F_advForce[k]/1000 << " "; // we need the output in kN
                morisonForce << data->eta[i*data->nx*data->ny+nn];
                morisonForce << "\n";


            }

        } else {

            morisonForce << data->t[i] << " ";

            Fout[i] = 0;

            double upper = 0;
            double lower = 0;
            Fout[i] = F[0]*(z_ocw[1]-z_ocw[0])/2;
            for (int k = 1;(unsigned)k<F.size()-1;k++){
                upper = (z_ocw[k+1]-z_ocw[k])/2;
                lower = (z_ocw[k]-z_ocw[k-1])/2;
                Fout[i] += F[k]*(upper+lower);
            }
            Fout[i] += F[F.size()]*(z_ocw[F.size()]-z_ocw[F.size()-1])/2;

            morisonForce << Fout[i]<< " "; // we need the output in kN
            morisonForce << data->eta[i*data->nx*data->ny+nn];
            morisonForce << "\n";



        }


    }

    morisonForce.close();

}

void advForce::interp_linear(const std::vector<double> f_ind,const std::vector<double> x_ind,const std::vector<double> x_out, std::vector<double>& f_out, int& id){

    double df;
    double dx;
    for (int i = 0;(unsigned)i<x_out.size();i++){

        double x = x_out[i];

        for (int j=0;(unsigned)j<x_ind.size()-1;j++){

            if (x<x_ind[0]){
                f_out[i] = f_out[0];
                break;
            }

            if (x_ind[j+1]>x){

                df = (f_ind[j+1]-f_ind[j]);
                dx = (x_ind[j+1]-x_ind[j]);

                f_out[i] = f_ind[j]+df/dx*(x-x_ind[j]);
                id =i;
                break;

            }
        }
    }

}

void advForce::readTable(){

    z_table.clear();z_table.resize(forceCoeffsAndDiameter->rowCount());
    D_table.clear();D_table.resize(forceCoeffsAndDiameter->rowCount());
    Cm_table.clear();Cm_table.resize(forceCoeffsAndDiameter->rowCount());
    Cd_table.clear();Cd_table.resize(forceCoeffsAndDiameter->rowCount());

    QDoubleSpinBox *sp;

    int zCol        = 0;
    int diameterCol = 1;
    int cdCol       = 2;
    int cmCol       = 3;

    for (int i =0;i<forceCoeffsAndDiameter->rowCount();i++){
        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,zCol);
//        z_table.insert(z_table.begin()+i,sp->value());
        z_table.at(i)= sp->value();

        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,diameterCol);
        D_table.at(i) =sp->value();

        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,cdCol);
        Cd_table.insert(Cd_table.begin()+i,sp->value());

        sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,cmCol);
        Cm_table.insert(Cm_table.begin()+i,sp->value());

    }

}

void advForce::plot_(){

    computeForce(xLocation->currentIndex());

    QVector<double> QV_t;
    QV_t.resize(data->nt);

    for (int i=1;i<data->nt-1;i++){
        QV_t[i] = data->t[i];
    }


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
    cPlot->graph()->setData(QV_t,Fout);
    cPlot->xAxis->setLabel("Time, t s");
    cPlot->yAxis->setLabel("Inline force, F N");
    QVector<double>::iterator Xmin = std::min_element(QV_t.begin(), QV_t.end());
    QVector<double>::iterator Xmax = std::max_element(QV_t.begin(), QV_t.end());
    QVector<double>::iterator Ymin = std::min_element(Fout.begin(), Fout.end());
    QVector<double>::iterator Ymax = std::max_element(Fout.begin(), Fout.end());
    cPlot->xAxis->setRange(*Xmin,*Xmax);
    cPlot->yAxis->setRange(*Ymin,*Ymax);
    cPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    cPlot->replot();

}

void advForce::save_(){

    if (forceCoeffsAndDiameter->rowCount()<2){
        QMessageBox msgBox;
        msgBox.setText("Please define a geometry");
        msgBox.exec();
    return;
    }

    if(advForceFile->text().length()==0){
        QMessageBox msgBox;
        msgBox.setText("Save using default name: \"OceanWave3D_advForce.wf\"");
        msgBox.setInformativeText("Do you want to continue?");
        msgBox.setStandardButtons(QMessageBox::Ok  | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Ok);

        int ret = msgBox.exec();

        switch (ret) {
        case QMessageBox::Ok:
            advForceFile->setText("OceanWave3D_advForce.wf");
            advForceFile_fileInfo.setFile("OceanWave3D_advForce.wf");

            computeForce(xLocation->currentIndex());


        case QMessageBox::Cancel:
            return;

        }
    }


    computeForce(xLocation->currentIndex());


}

void advForce::saveGeometry(){

    if(geometryFile->text().length()==0){
        QMessageBox msgBox;
        msgBox.setText("Save using default name: \"pile.monopile\"");
        msgBox.setInformativeText("Do you want to continue?");
        msgBox.setStandardButtons(QMessageBox::Ok  | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Ok);

        int ret = msgBox.exec();

        switch (ret) {
        case QMessageBox::Ok:
            geometryFile->setText("pile.monopile");
            geometryFile_fileInfo.setFile("pile.monopile");

        case QMessageBox::Cancel:
            return;

        }
    }

    QDoubleSpinBox *sp;
    QFile geoFile(geometryFile_fileInfo.fileName());
    geoFile.open(QIODevice::WriteOnly);
    QTextStream stream(&geoFile);
    // write header

    stream << geometryHeader->text();
    if (!geometryHeader->text().endsWith("\n")){stream << "\n";}

    // write number of rows
    stream << forceCoeffsAndDiameter->rowCount() << "\n";

    // write the data
    for (int i =0;i<forceCoeffsAndDiameter->rowCount();i++){
        for (int j=0;j<forceCoeffsAndDiameter->columnCount();j++){
            sp = (QDoubleSpinBox*)forceCoeffsAndDiameter->cellWidget(i,j);
            stream << sp->value() << " ";
        }

        stream << "\n";


    }
    geoFile.close();
}
