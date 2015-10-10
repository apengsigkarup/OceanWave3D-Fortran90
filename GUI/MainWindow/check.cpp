#include "mainwindow.h"
#include "checkdialog.h"
#include "ui_mainwindow.h"
#include <QMessageBox>

bool MainWindow::checkCase()
{
    bool finalCheck = true;
    checkDialog *check = new checkDialog(this);
    double L;
    // Spatial resolution
    //
bool irr;
    // wave length
    if (ui->waveType->currentIndex()==1){ // stream function wave

        irr = false;

        if (ui->LorP_ComboBox->currentIndex()==0){ // Wave length specified
            L = ui->SF_T->value();
        }
        else { // wave period
            if (ui->SF_h->value()==0){errorMsgSFhZero();return false;}
            double k = dispersion_T(ui->SF_h->value(),ui->SF_T->value());
            L = 2*M_PI/k;
        }


    } else if ((ui->waveType->currentIndex()==2) | (ui->waveType->currentIndex()==3)) { // JONSWAP & PM

        irr = true;

        L = 2*3.1415/ui->maxkh->value()*ui->h->value();

    }
    else {
        errorMsgParameterCheck();
        return false;

    }

    // dx
    double dx = ui->length->value()/(ui->nx->value()-1);

    // ppwl
    double spatialResolution = L/dx+1;

    check->spatialResolution(spatialResolution,irr) ? finalCheck=true  : finalCheck = false;

    // Temporal resolution
    //
    double T;
    if (ui->waveType->currentIndex()==1){ // stream function wave

        if (ui->LorP_ComboBox->currentIndex()==1){ // Wave period
            T = ui->SF_T->value();
        }
        else { // wave length
            if (ui->SF_h->value()==0){errorMsgSFhZero();return false;}
            double omega = dispersion_L(ui->SF_h->value(),ui->SF_T->value());
            T = 2*M_PI/omega;
        }


    } else if ((ui->waveType->currentIndex()==2) | (ui->waveType->currentIndex()==3)) { // JONSWAP & PM
        L = 2*M_PI/ui->maxkh->value()*ui->depth->value();
        double omega = dispersion_L(ui->h->value(),L);
        T = 2*M_PI/omega;


    }
    else {
        errorMsgParameterCheck();
        return false;

    }

    check->temporalResolution(T/ui->dt->value(),irr) ? finalCheck=true : finalCheck=false;


    // We don't want to base the relaxation zones on the sortest waves
    // in the irregular spectrum - so we recompute based on the peak period
    if ((ui->waveType->currentIndex()==2) | (ui->waveType->currentIndex()==3)){
        L = 2*M_PI/dispersion_T(ui->h->value(),ui->Tp->value());
    }
    // Generation Zone
    //
    check->generationZone(ui->xGenStart->value(),ui->xGenEnd->value(),L) ? finalCheck=true : finalCheck=false;

    // Generation Zone
    //
    check->relaxationZone(ui->xAbsorbStart->value(),ui->xAbsorbEnd->value(),L,ui->length->value()) ? finalCheck=true  : finalCheck = false;

    // depth
    //
    double depth;
    irr ? depth = ui->h->value() : depth = ui->SF_h->value();

    check->depth(ui->depth->value(),depth) ? finalCheck=true  : finalCheck == false;

    check->show();
    if (finalCheck == false){
        return false;
    } else {
        return true;
    }
}
