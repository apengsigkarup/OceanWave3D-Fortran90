#include "externaloutput.h"
#include "ui_externaloutput.h"

externalOutput::externalOutput(const convert &data_, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::externalOutput)
{
    ui->setupUi(this);
}

externalOutput::~externalOutput()
{
    delete ui;
}
