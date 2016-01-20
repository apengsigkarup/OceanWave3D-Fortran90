#ifndef EXTERNALOUTPUT_H
#define EXTERNALOUTPUT_H

#include <QWidget>
#include "convert.h"

namespace Ui {
class externalOutput;
}

class externalOutput : public QWidget
{
    Q_OBJECT
    
public:
    explicit externalOutput(const convert& data_,QWidget *parent = 0);
    ~externalOutput();
    
private:
    Ui::externalOutput *ui;
};

#endif // EXTERNALOUTPUT_H
