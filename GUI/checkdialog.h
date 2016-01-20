#ifndef CHECKDIALOG_H
#define CHECKDIALOG_H

#include <QDialog>
#include "math.h"
#include "versions.h"

namespace Ui {
class checkDialog;
}



class checkDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit checkDialog(QWidget *parent = 0);
    ~checkDialog();
    bool spatialResolution(double resolution, bool irr);
    bool temporalResolution(double resolution, bool);
    bool generationZone(double, double, double );
    bool relaxationZone(double, double, double, double );
    bool depth(const double, const double);

private:
    Ui::checkDialog *ui;

};




#endif // CHECKDIALOG_H
