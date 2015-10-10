#include "checkdialog.h"
#include "QPainter"
#include "ui_checkdialog.h"

checkDialog::checkDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::checkDialog)
{
    ui->setupUi(this);

}

checkDialog::~checkDialog()
{
    delete ui;
}



bool checkDialog::spatialResolution(double resolution, bool irr)
{

    bool finalCheck = true;
    QPixmap pixmapGood(":/pictures/checkMark.png");
    QPixmap pixmapBad(":/pictures/bad.png");
    QPixmap pixmapMed(":/pictures/warning.png");

    QString outstring = "Points per wave length = " + QString::number(resolution,'f',2);
    ui->spatialRes_label->setText(outstring);

    if (!irr){
    if (resolution <=5)
    {
        ui->spatialRes_image->setPixmap(pixmapBad);
        ui->explanation_spatial->setText("Insufficient points per wave length");
        finalCheck = false;
    } else if ((resolution>5) & (resolution<15))
    {
        ui->spatialRes_image->setPixmap(pixmapMed);
        ui->explanation_spatial->setText("limited points per wave length");
        finalCheck = false;
    } else {
        ui->spatialRes_image->setPixmap(pixmapGood);
    }
    } else {
        if (resolution <5)
        {
            ui->spatialRes_image->setPixmap(pixmapBad);
            ui->explanation_spatial->setText("Insufficient points per wave length");
            finalCheck = false;
        } else {
            ui->spatialRes_image->setPixmap(pixmapGood);
          ui->explanation_spatial->setText("Note: For irregular waves we are less\n restrictive");
        }

    }
    ui->spatialRes_image->setScaledContents(true);

    return finalCheck;
}


bool checkDialog::temporalResolution(double resolution,bool irr)
{



    bool finalCheck = true;
    QPixmap pixmapGood(":/pictures/checkMark.png");
    QPixmap pixmapBad(":/pictures/bad.png");
    QPixmap pixmapMed(":/pictures/warning.png");

    QString outstring = "Points per wave period = " + QString::number(resolution,'f',2);
    ui->temporalRes_label->setText(outstring);
    if (!irr){
        if (resolution <=20)
        {
            ui->temporalRes_image->setPixmap(pixmapBad);
            ui->explanation_temporal->setText("Insufficient points per wave period");
            finalCheck = false;
        } else if ((resolution>20) & (resolution<30))
        {
            ui->temporalRes_image->setPixmap(pixmapMed);
            ui->explanation_temporal->setText("Limited points per wave period");
            finalCheck = false;
        } else {
            ui->temporalRes_image->setPixmap(pixmapGood);
        }
    } else {
        if (resolution <=20)
        {
            ui->temporalRes_image->setPixmap(pixmapBad);
            ui->explanation_temporal->setText("Insufficient points per wave period");
            finalCheck = false;
        } else {
            ui->temporalRes_image->setPixmap(pixmapGood);
            ui->explanation_temporal->setText("Note: For irregular waves we are less\n restrictive");
        }


    }
    ui->temporalRes_image->setScaledContents(true);

    return finalCheck;
}

bool checkDialog::generationZone(double start, double end, double L){

    bool finalCheck = true;
    double length = (end - start)/L;

    QPixmap pixmapGood(":/pictures/checkMark.png");
    QPixmap pixmapBad(":/pictures/bad.png");
    QString outstring;
    QString explanation = "";

    outstring = "start = " + QString::number(start,'f',2);
    ui->genZone_start_label->setText(outstring);
    bool startCheck = start ==0 ? true : false;
    if (!startCheck){explanation="Check the start point. \n";}

    outstring = "end = " + QString::number(end,'f',2);
    ui->genZone_end_label->setText(outstring);
    bool endCheck = end>start ? true : false;
    if (!endCheck){explanation+="Check the end point. \n";}

    outstring = "length = " + QString::number(length,'f',2) + " (normalized by L)";
    ui->genZone_length_label->setText(outstring);

    bool lengthCheck;
    if ((length<0.5) | (length > 5)){
        lengthCheck = false;
        explanation += "The length is incorrect \n";
    }  else {
        lengthCheck = true;
    }

    if (!startCheck | !endCheck | !lengthCheck )
    {
        ui->genZone_image->setPixmap(pixmapBad);
        finalCheck = false;
    } else {
        ui->genZone_image->setPixmap(pixmapGood);
    }

     ui->genZone_image->setScaledContents(true);
     ui->explanation_genZone->setText(explanation);

     return finalCheck;
}


bool checkDialog::relaxationZone(double start, double end, double L, double domainLength){
    bool finalCheck = true;
    double length = (end - start)/L;

    QPixmap pixmapGood(":/pictures/checkMark.png");
    QPixmap pixmapBad(":/pictures/bad.png");
    QString outstring;
    QString explanation = "";

    outstring = "start = " + QString::number(start,'f',2);
    ui->endZone_start_label->setText(outstring);
    bool startCheck = start < end ? true : false;
    if (!startCheck){explanation="Check the start point. \n";}

    outstring = "end = " + QString::number(end,'f',2);
    ui->endZone_end_label->setText(outstring);
    bool endCheck = end==domainLength ? true : false;
    if (!endCheck){explanation+="Check the end point. \n";}

    outstring = "length = " + QString::number(length,'f',2) + " (normalized by L)";
    ui->endZone_length_label->setText(outstring);

    bool lengthCheck;

    if ((length<0.9) | (length > 5)){
        lengthCheck = false;
        explanation += "The length is incorrect \n";
    } else {
        lengthCheck = true;
    }

    if (!startCheck | !endCheck | !lengthCheck)
    {
        ui->endZone_image->setPixmap(pixmapBad);
        finalCheck = false;

    } else {
        ui->endZone_image->setPixmap(pixmapGood);
    }
     ui->endZone_image->setScaledContents(true);
     ui->explanation_endZone->setText(explanation);

     return finalCheck;
}

bool checkDialog::depth(const double h_geo, const double h_waveT){
    bool finalCheck;
    QPixmap pixmapGood(":/pictures/checkMark.png");
    QPixmap pixmapWarning(":/pictures/warning.png");

    ui->depth_geo->setText("Geometry: depth = " + QString::number(h_geo,'f',2) + " m");
    ui->depth_waveTheory->setText("Wave theory: depth = " + QString::number(h_waveT,'f',2) + " m");

    fabs(h_geo-h_waveT)<0.00001 ? finalCheck = true : finalCheck = false;

    if (finalCheck){
        ui->depth_image->setPixmap(pixmapGood);
    } else {
        ui->depth_image->setPixmap(pixmapWarning);
    }
    ui->depth_image->setScaledContents(true);

    return finalCheck;
}
