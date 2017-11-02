#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "math.h"
#include <cmath>
#include <QProcess>


void MainWindow::run(){

    //    QTextStream out(stdout);
    //     out << QString("Some text");
    writeInputFile();
    QProcess *OCW = new QProcess();
    OCW->start("/usr/bin/xterm -hold -e \"OceanWave3D\"");
}

void MainWindow::openWorkDirDialog(){

    dialog.setFileMode(QFileDialog::Directory);
    dialog.setOption(QFileDialog::ShowDirsOnly);
    dialog.exec();

    dir.setCurrent(dialog.directory().path());

    ui->workingDir->setText(dialog.directory().path());
}


void MainWindow::clearCase(){

    QString argument("/bin/rm -f fort.* Kinematics*");
    std::system(qPrintable(argument));


}
void MainWindow::openFile(){
    QString fileName = QFileDialog::getOpenFileName(
                this,
                "Select file",
                dir.currentPath(),
                "OCW3D (*.inp);;"
                );


    QFile inputFile(fileName);

    QString tmp_line;
    QStringList tmp_list;
    QString errMsg;

    if (!inputFile.open(QIODevice::ReadOnly)){
        errorMsgUnknownFile();
        return;
    }

    QFileInfo fileInfo(inputFile.fileName());


    // set working directory
    //
    dir.setCurrent(fileInfo.path());
    ui->workingDir->setText(fileInfo.path());

    // Header line
    //
    tmp_line = inputFile.readLine();tmp_line.remove(QRegExp("[\\n\\t\\r]"));
    ui->header_input->setText(tmp_line);

    // Second line: initial conditions, wave type, acceleration.
    //
    readAndSplit(inputFile,tmp_list);


    //    if (tmp_list.size()>=1){ui->initialCondition->setCurrentIndex(tmp_list[0].toInt());}
    if (tmp_list.size()>=2){
        if (tmp_list[1].toInt()==1) {ui->waveType->setCurrentIndex(1);}
        if (tmp_list[1].toInt()==2) {ui->waveType->setCurrentIndex(2);} // we update with JONSWAP or PM later
        if (tmp_list[1].toInt()==3) {ui->waveType->setCurrentIndex(5);}
    }
    if (tmp_list.size()>=3){
        if ((tmp_list[2].toDouble()>100) | (tmp_list[2].toDouble()==0)){
            tmp_list[2] = "0";
            ui->breakingWidget->setEnabled(false);
            ui->checkBox_breakingWidget->setChecked(true);
        } else {
            ui->breakingWidget->setEnabled(true);
            ui->checkBox_breakingWidget->setChecked(false);
        }
        ui->breaking_beta0->setValue(tmp_list[2].toDouble());

    }

    // third line
    //
   readAndSplit(inputFile,tmp_list);

    ui->length->setValue(tmp_list[0].toDouble());
    ui->width->setValue(tmp_list[1].toDouble());
    ui->depth->setValue(tmp_list[2].toDouble());
    ui->nx->setValue(tmp_list[3].toInt());
    ui->ny->setValue(tmp_list[4].toInt());
    ui->nz->setValue(tmp_list[5].toInt());

    double dx = tmp_list[0].toDouble()/(tmp_list[3].toDouble()-1);
    double dy;
    ui->ny->value()==1 ? dy = tmp_list[1].toDouble() : dy = tmp_list[1].toDouble()/(tmp_list[4].toDouble()-1);
    ui->geometryType->setCurrentIndex(0);
    tmp_list[8].toInt()>0 ? ui->sz->setChecked(true):ui->sz->setChecked(false);

    if (tmp_list.length()>12){
        ui->selectedGridFile->setText(tmp_list[12].toLatin1());
        ui->depth->setValue(-1*tmp_list[2].toDouble());
        error_grid_checkDomainLength(); // throw warning to user

        ui->geometryType->setCurrentIndex(1);
        ui->output->setCurrentIndex(1);
        ui->length->setEnabled(true);
    }


    // forth line, alpha, beta, gamma, a b c
    //
    readAndSplit(inputFile,tmp_list);

    ui->alpha->setValue(tmp_list[0].toInt());
    ui->beta->setValue(tmp_list[1].toInt());
    ui->gamma->setValue(tmp_list[2].toInt());
    ui->a->setValue(tmp_list[3].toInt());
    ui->b->setValue(tmp_list[4].toInt());
    ui->c->setValue(tmp_list[5].toInt());


    // fith line: nt dt
    //
    readAndSplit(inputFile,tmp_list);
    double time = tmp_list[0].toDouble()*tmp_list[1].toDouble();
    double dt = tmp_list[1].toDouble(); // we need this value later

    ui->timeDuration->setValue(time);
    ui->dt->setValue(tmp_list[1].toDouble());

    // sixth line: g \rho
    //
    readAndSplit(inputFile,tmp_list);
    if (tmp_list[0].toDouble()==0){
        ui->gravity_input->setValue(9.81);
    } else {
        ui->gravity_input->setValue(tmp_list[0].toDouble());
    }

    if (tmp_list[1].toDouble()==0) {
        ui->Density->setValue(1025);
    } else {
        ui->Density->setValue(tmp_list[1].toDouble());
    }

    // seventh line: numerical settings (SKIP)
    //
    tmp_line = inputFile.readLine();

    // stream function solution
    //
    readAndSplit(inputFile,tmp_list);
    if (ui->waveType->currentIndex()==1){
        ui->SF_H->setValue(tmp_list[0].toDouble());
        ui->SF_h->setValue(tmp_list[1].toDouble());
        if (tmp_list[4].toInt()==0){
            ui->LorP_ComboBox->setCurrentIndex(0);
            ui->SF_TLabel->setText("Wave length");
            ui->SF_TL_unit->setText("m");
            ui->SF_T->setValue(tmp_list[2].toDouble());
        } else {
            ui->LorP_ComboBox->setCurrentIndex(1);
            ui->SF_TLabel->setText("Wave period");
            ui->SF_TL_unit->setText("s");
            ui->SF_T->setValue(tmp_list[3].toDouble());
        }
        ui->SF_U->setValue(tmp_list[5].toDouble());
        ui->stokesOrEuler->setCurrentIndex(tmp_list[6].toInt());
        ui->SF_n->setValue(tmp_list[8].toInt());
    }

    tmp_line.clear();
    tmp_list.clear();

    // outputs
    readAndSplit(inputFile,tmp_list);

    if ((tmp_list[0].toInt()==0) & (tmp_list[1].toInt()==0)){ // no outputs
        ui->storeAscii_onOff->setChecked(false);
        ui->nOutFiles->setValue(0);
    } else {
        if (tmp_list[0].toInt()<0){// we store ascii files
            ui->storeAscii_onOff->setChecked(true);
            ui->ACCII_label->setVisible(true);
            ui->ASCII_label2->setVisible(true);
            ui->nASCII->setVisible(true);
            ui->nASCII->setValue(-1*tmp_list[0].toInt());
        } else {
            ui->storeAscii_onOff->setChecked(false);
            ui->ACCII_label->setVisible(false);
            ui->ASCII_label2->setVisible(false);
            ui->nASCII->setVisible(false);

        }

        if (tmp_list[3].toInt()>0){// Kinematics or wave gauge files

            int nFiles = tmp_list[3].toInt();
            on_outputWidgetChanged(nFiles);


            if (tmp_list[1].toInt()==20 ){
                ui->DropDownListOutputType->setCurrentIndex(1);
            } else {
                ui->DropDownListOutputType->setCurrentIndex(2);
            }

            ui->tableWidget->setVisible(true);
            ui->nOutFiles->setValue(nFiles);
            ui->tableWidget->setRowCount(0);


            QDoubleSpinBox *sp;


            for (int i =0; i<nFiles;i++){
                readAndSplit(inputFile,tmp_list);

                tmp_list.removeAt(2);
                tmp_list.removeAt(5);
                QList<double> outPutList;



                ui->tableWidget->insertRow(i);

                if (ui->DropDownListOutputType->currentIndex()==1) {
                    outPutList.push_back((tmp_list[0].toDouble()-1)*dx);
                    outPutList.push_back((tmp_list[1].toDouble()-1)*dx);
                    outPutList.push_back((tmp_list[2].toDouble()-1)*dy);
                    outPutList.push_back((tmp_list[3].toDouble()-1)*dy);
                    outPutList.push_back((tmp_list[4].toDouble()-1)*dt);
                    outPutList.push_back((tmp_list[5].toDouble())*dt);


                    for (int j = 0; j < 6 ; j++) {
                        ui->tableWidget->setCellWidget(i,j,new QDoubleSpinBox(ui->tableWidget));
                        sp = (QDoubleSpinBox*)ui->tableWidget->cellWidget(i,j);
                        sp->setMaximum(99999999);
                        sp->setValue(outPutList[j]);
                    }
                }
                if (ui->DropDownListOutputType->currentIndex()==2){

                    outPutList.push_back((tmp_list[0].toDouble()-1)*dx);
                    outPutList.push_back((tmp_list[2].toDouble()-1)*dy);
                    outPutList.push_back((tmp_list[4].toDouble()-1)*dt);
                    outPutList.push_back((tmp_list[5].toDouble())*dt);

                    if (ui->ny->value()==1) {outPutList[1] =1; }
                    for (int j = 0; j < 2 ; j++) {
                        ui->tableWidget->setCellWidget(i,j,new QDoubleSpinBox(ui->tableWidget));
                        sp = (QDoubleSpinBox*)ui->tableWidget->cellWidget(i,j);
                        sp->setMaximum(99999999);
                        sp->setValue(outPutList[j]);
                    }

                }


            }

        }

    }


    // linear or nonlinear
    //
    readAndSplit(inputFile,tmp_list);

    if (tmp_list[0].toInt()==1){
        ui->nonlin_onOff->setCurrentIndex(0);
    } else {
        ui->nonlin_onOff->setCurrentIndex(1);

    }

    // boudary smoothing (SKIP)
    //
    tmp_line = inputFile.readLine();

    // relaxation zones
    //
    readAndSplit(inputFile,tmp_list);

    if (tmp_list[0].toInt()==1){

        ui->rampTime->setValue(tmp_list[1].toDouble());
        int relaxZone(tmp_list[2].toInt());

        // we read the first zone
        readAndSplit(inputFile,tmp_list);

        ui->xGenStart->setValue(tmp_list[0].toDouble());
        ui->xGenEnd->setValue(tmp_list[1].toDouble());
        ui->yGenStart->setValue(tmp_list[2].toDouble());
        ui->yGenEnd->setValue(tmp_list[3].toDouble());

        if (relaxZone==2){
             readAndSplit(inputFile,tmp_list);
            ui->xAbsorbStart->setValue(tmp_list[0].toDouble());
            ui->xAbsorbEnd->setValue(tmp_list[1].toDouble());
            ui->yAbsorbStart->setValue(tmp_list[2].toDouble());
            ui->yAbsorbEnd->setValue(tmp_list[3].toDouble());
        }
    }


    // pressure damping zones
    //
    readAndSplit(inputFile,tmp_list);

    if (tmp_list[0].toInt()==1){// pressure damping is on
        ui->pressureDampingOrRelax->setCurrentIndex(0);
        tmp_line = inputFile.readLine();tmp_line.remove(QRegExp("[\\n\\t\\r]"));
        tmp_list = tmp_line.split(" ");

        ui->xAbsorbStart->setValue(tmp_list[0].toDouble());
        ui->xAbsorbEnd->setValue(tmp_list[1].toDouble());
        ui->yAbsorbStart->setValue(tmp_list[2].toDouble());
        ui->yAbsorbEnd->setValue(tmp_list[3].toDouble());

    } else {
        ui->pressureDampingOrRelax->setCurrentIndex(1);
    }


    // SWENSE (SKIP)
    //
    tmp_line = inputFile.readLine();

    // curvelinear (SKIP)
    //
    tmp_line = inputFile.readLine();

    // irregular wave parameters
    //
    if (ui->waveType->currentIndex()==2){

        readAndSplit(inputFile,tmp_list);

        if (tmp_list[0].toInt() == 0) {
            ui->waveType->setCurrentIndex(3); // PM
        } else if (tmp_list[0].toInt() == 1) {
            ui->waveType->setCurrentIndex(2); // jonswap old syntax
        } else if (tmp_list[0].toInt() == 2) {
            ui->waveType->setCurrentIndex(4);
        } else if (tmp_list[0].toInt() == 3) {
            ui->waveType->setCurrentIndex(2); // jonswap new syntax
        }
    }
    if ((ui->waveType->currentIndex()==2) | (ui->waveType->currentIndex()==3) )
    {
        ui->Hs->setValue(tmp_list[2].toDouble());
        ui->Tp->setValue(tmp_list[1].toDouble());
        ui->h->setValue(tmp_list[3].toDouble());
        ui->maxkh->setValue(tmp_list[4].toDouble());
        ui->seed->setValue(-1*tmp_list[5].toInt());

        if (tmp_list[0].toInt() == 1){
            ui->gamma_jonswap->setValue(3.3); // JONSWAP old syntax
        } else {
            ui->gamma_jonswap->setValue(tmp_list[9].toDouble()); // JONSWAP with variable gamma
        }



    }
    if (ui->waveType->currentIndex()==4) {
        ui->Hs->setValue(tmp_list[2].toDouble());
        ui->Tp->setValue(tmp_list[1].toDouble());
        ui->h->setValue(tmp_list[3].toDouble());
        ui->maxkh->setValue(tmp_list[4].toDouble());
        ui->seed->setValue(-1*tmp_list[5].toInt());
        ui->irr_x0->setValue(tmp_list[7].toDouble());
        ui->irr_y0->setValue(tmp_list[8].toDouble());
        QFileInfo waveFile(tmp_list[9].toLatin1());
        ui->selectedWaveFile_eta->setText(waveFile.fileName());
    }
    if (ui->waveType->currentIndex()==5){
        tmp_line = inputFile.readLine();tmp_line.remove(QRegExp("[\\n\\t\\r]"));
        tmp_list = tmp_line.split(" ");
        QFileInfo waveFile(tmp_list[2].toLatin1());
        ui->rampTime_waveFile->setValue(tmp_list[0].toDouble());
        ui->selectedWaveFile->setText(waveFile.fileName());
    }

}



void MainWindow::writeInputFile()
{
    // Read inputs from tab GENERAL
    // Read header
    QString header = ui->header_input->text();
    // Read mode
    int mode;
    ui->nonlin_onOff->currentIndex()==0 ? mode=1 : mode = 0;

    // Gravity
    double gravity = ui->gravity_input->value();
    double Density = ui->Density->value();
    // Breaking filter
    double breaking = ui->breaking_beta0->value();
    if (breaking==0){breaking = 1000;}
    // Time duration
    double timeDur = ui->timeDuration->value();
    tStart = 0; // ui->tStart->value(); We have removed this feature since it adds more confusion than functionality

    // Read initial conditions
    int initialConditions = 0; // still water is hardcoded. Alternatives are 0) still water, 1) standing wave

    // read discretization
    int alpha = ui->alpha->value();
    int beta = ui->beta->value();
    int gamma = ui->gamma->value();
    int a = ui->a->value();
    int b = ui->b->value();
    int c = ui->c->value();
    double dt = ui->dt->value();

    // Compute the number of time steps
    int Nsteps = timeDur/dt;


    // Read geometry
    double length = ui->length->value();
    double width = ui->width->value();
    double depth = ui->depth->value();
    int nx = ui->nx->value();
    int ny = ui->ny->value();
    int nz = ui->nz->value();
    bool sx = false; // This option is not supported by the GUI
    bool sy = false; // This option is not supported by the GUI
    bool sz = ui->sz->isChecked();
    double epsilon = 1e-12;
    double dx = length/(nx-1);
    double dy = (width+epsilon)/(ny-1);

    // Read wave generation
    // wave theory
    int waveTheory;

    if (ui->waveType->currentIndex()==1){ // Stream Fucntion
        SF_H = ui->SF_H->value();
        SF_T = ui->SF_T->value();
        SF_h = ui->SF_h->value();
        SF_U = ui->SF_U->value();
        SorE = ui->stokesOrEuler->currentIndex();
        SF_n = ui->SF_n->value();
        LorT = ui->LorP_ComboBox->currentIndex();
        SF_L=100;//Dummy value
        waveTheory = 1;
        if (LorT==0) {SF_L =SF_T;}
    }
    if ((ui->waveType->currentIndex()==2) | (ui->waveType->currentIndex()==3)){ //JONSWAP or PM
        irrFilename = QString(' ');
        Hs = ui->Hs->value();
        Tp = ui->Tp->value();
        h = ui->h->value();
        gamma_jonswap = ui->gamma_jonswap->value();
        seed = ui->seed->value();
        khmax = ui->maxkh->value();
        waveTheory=2;
        ui->waveType->currentIndex()==2 ? i_spec=3 : i_spec=0;
    }
    if (ui->waveType->currentIndex()==5){
        waveTheory=3;
    }

    if (ui->waveType->currentIndex()==4) {
        irrFilename = ui->selectedWaveFile_eta->text();
        Hs = ui->Hs->value();
        Tp = ui->Tp->value();
        h = ui->h->value();
        //        gamma = ui->gamma->value();
        seed = ui->seed->value();
        khmax = 2*M_PI/dx*h;
        waveTheory=2;
        i_spec = 2;
    }

    // Read relaxation / damping zones
    //Generation
    double xGen[2];double yGen[2];
    xGen[0] = ui->xGenStart->value();xGen[1] = ui->xGenEnd->value();
    yGen[0] = ui->yGenStart->value();yGen[1] = ui->yGenEnd->value();
    double rampTime = ui->rampTime->value();

    double xAbsorb[2];double yAbsorb[2];
    xAbsorb[0] = ui->xAbsorbStart->value();xAbsorb[1] = ui->xAbsorbEnd->value();
    yAbsorb[0] = ui->yAbsorbStart->value();yAbsorb[1] = ui->yAbsorbEnd->value();


    // Read output files
    bool storeAscii = ui->storeAscii_onOff->isChecked();
    int outputType = ui->DropDownListOutputType->currentIndex();
    int nCol;
    std::vector<double> resolution;
    if (outputType==1) {
        nCol = 6;
        double Res[] ={dx,dx,dy,dy,dt,dt};
        resolution.insert(resolution.end(), &Res[0], &Res[sizeof(Res) / sizeof(double)]);

    } else if (outputType==2) {
        nCol = 2;
        double Res[] ={dx,dy};
        resolution.insert(resolution.end(), &Res[0], &Res[sizeof(Res) / sizeof(double)]);
    }
    int nASCIIsteps = ui->nASCII->value();


    QDoubleSpinBox* sp;
    Double2d outputValues(boost::extents[ui->tableWidget->rowCount()][6]);
    double tmp;
    for (int i=0;i<ui->tableWidget->rowCount();i++){
        for (int j = 0; j < nCol ; j++) {
            sp = (QDoubleSpinBox*)ui->tableWidget->cellWidget(i,j);
            tmp = round(sp->value()/resolution[j]+1);
            outputValues[i][j]= tmp;
        }

        if (outputType==2) {
            outputValues[i][5] = Nsteps;
            outputValues[i][4] = 1;
            outputValues[i][3] = outputValues[i][1];
            outputValues[i][2] = outputValues[i][1];
            outputValues[i][1] = outputValues[i][0];
        }

        if (outputValues[i][4]==0){outputValues[i][4]=1;}
        if (outputValues[i][5]>Nsteps){outputValues[i][5]=Nsteps;}
        if (ny==1) {outputValues[i][2]=1;outputValues[i][3]=1;}

    }

    //     Write OceanWave3D input file
    QFile myfile("OceanWave3D.inp");
    myfile.open (QIODevice::WriteOnly);
    QTextStream outStream(&myfile);
    outStream << header;
    if (!header.endsWith("\n")){outStream << "\n";}

    outStream << initialConditions << " " << waveTheory << " " << breaking << "\n";

    if (ui->geometryType->currentIndex()==0) {
        outStream << length << " "<< width << " " << depth << " " << nx << " " << ny << " " << nz << " " << sx << " " << sy << " " << sz << " 1 1 1\n";
    } else {
        QString gridName = ui->selectedGridFile->text();
        outStream << length << " "<< width << " " << -1*depth << " " << nx << " " << ny << " " << nz << " " << sx << " " << sy << " " << sz << " 1 1 1 " << gridName;
        if (!gridName.endsWith("\n")){outStream << "\n";}
    }
    outStream << alpha << " " << beta << " " << gamma << " " << a << " " << b << " " << c << "\n";
    outStream << Nsteps << " " << dt << " " << "1 0 0 " << tStart << "\n";
    outStream << gravity <<" " << Density << "\n";
    outStream << "1 1 0 23 1e-8 1e-6 1 V 1 1 2 \n";
    if (waveTheory==1){
        outStream << SF_H << " " << SF_h << " " << SF_L << " " << SF_T << " "<< LorT << " " << SF_U << " " << SorE << " " << "8 " << SF_n << "\n";
    } else { outStream << "1 1 1 1 1 1 1 1 1\n";}

    if ((ui->DropDownListOutputType->currentIndex()==0)&(~storeAscii)){ // No outputs requested
        outStream << "0 0 \n";
    } else { // Outputs are requested

        if ((outputType==1)&storeAscii){ // Kinematics and ascii files
            outStream << -1*nASCIIsteps << " 20 0 " << ui->nOutFiles->value() << "\n";
        }
        if ((outputType==2)&storeAscii){ // Wave gagues and ascii files
            outStream << -1*nASCIIsteps << " 30 0 " << ui->nOutFiles->value() << "\n";
        }
        if ((outputType==0)&(storeAscii)){ // Only ASCII files
            outStream << -1*nASCIIsteps << " 1 0 0\n";
        }
        if ((outputType==1)&(~storeAscii)){ // Only Kinematics files
            outStream << "0 20 1 " << ui->nOutFiles->value() << "\n";
        }
        if ((outputType==2)&(~storeAscii)){ // Only wave gagues
            outStream << "0 30 30 " << ui->nOutFiles->value() << "\n";
        }


        if (outputType>0){ // Writing output positions for binary files and wave guages

            for (int i = 0;i<ui->nOutFiles->value();i++){ // loop over rows
                for (int j = 0; j<6;j=j+2){
                    outStream << outputValues[i][j] << " " << outputValues[i][j+1] << " 1 ";
                }
                outStream << "\n";
            }

        }

    }

    outStream << mode << " 0\n"; // We hardcode the surface pressure term
    if (ui->waveType->currentIndex()==5) {
        outStream << "1 6 10 1 1 1 \n"; // for wave padle we need stronger boundary filtering
    } else {
        outStream << "1 6 10 0.08 0.08 0.4 \n"; // SG-filtering
    }


    if ((ui->waveType->currentIndex()==5) & (ui->pressureDampingOrRelax->currentIndex()==0)){
        outStream << "0 " << rampTime << " 0 X 0 \n"; // For wave paddle signal and pressure zone we don't need any relaxation zones
        outStream << "1 1\n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 1 1 0\n";
    } else if ((ui->waveType->currentIndex()==5) & (ui->pressureDampingOrRelax->currentIndex()==1)){ // Only absorbing relaxation zone
        outStream << "1 " << rampTime << " 1 X 0 \n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 9 3.5 X 0 X 0 \n";
        outStream << "0 0\n";
    } else if ((ui->waveType->currentIndex()!=5) & (ui->pressureDampingOrRelax->currentIndex()==0)) {
        outStream << "1 " << rampTime << " 1 X 0 \n";
        outStream << xGen[0] << " " << xGen[1] << " " << yGen[0] << " " << yGen[1] << " -9 3.5 X 1 X 0 \n";
        outStream << "1 1\n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 1 1 0\n";
    } else if ((ui->waveType->currentIndex()!=5) & (ui->pressureDampingOrRelax->currentIndex()==1)) {
        outStream << "1 " << rampTime << " 2 X 0 \n";
        outStream << xGen[0] << " " << xGen[1] << " " << yGen[0] << " " << yGen[1] << " -9 3.5 X 1 X 0 \n";
        outStream << xAbsorb[0] << " " << xAbsorb[1] << " " << yAbsorb[0] << " " << yAbsorb[1] << " 9 3.5 X 0 X 0 \n";
        outStream << "0 0\n";
    }

    outStream << "0 0 0 0 0 0 0\n"; // SWENSE
    outStream << "0 \n";
    // Irregualr waves`
    if (waveTheory==2){
        if (i_spec==2){ // Irregular wave with input file
            outStream << i_spec << " " << Tp << " " << Hs << " " << depth << " " << khmax << " " << -1*seed << " " << -1*seed << " " << ui->irr_x0->value() << " " << ui->irr_y0->value() << " " << irrFilename << " \n";
        }
        else if (i_spec==3) { // irregular wave with defined gamma vlaue
            outStream << i_spec << " " << Tp << " " << Hs << " " << depth << " " << khmax << " " << -1*seed << " " << -1*seed << " " << ui->irr_x0->value() << " " << ui->irr_y0->value() << " " << gamma_jonswap << " \n";
        }
    }
    if (waveTheory==3){ outStream << ui->rampTime_waveFile->value() << " " << 1 <<  " " << ui->selectedWaveFile->text(); }
    //
    myfile.close();
    ui->run->setEnabled(true);
}



void MainWindow::readAndSplit(QFile &argIn,QStringList &argOut){

    argOut.clear();

    QString line;
    line = argIn.readLine();
    line.remove(QRegExp("[\\n\\t\\r]"));
    argOut = line.split(" ");
    for (int i = 0;i<argOut.size();i++){
        if (argOut[i].isEmpty()){
            argOut.removeAt(i);
        }
    }

}
