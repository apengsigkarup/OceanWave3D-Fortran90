#include "mainwindow.h"
#include "ui_mainwindow.h"


// Select grid file
//

void MainWindow::selectGridFile()
{
QString fileName = QFileInfo(QFileDialog::getOpenFileName(
            this,
            "Select file",
            dir.currentPath(),
            "All files (*.*);;"
            )).fileName();

    ui->selectedGridFile->setText(fileName);

}
// Change geometry type; Simple gometry, flat bed;Custom geometry (2D)
//
void MainWindow::geometryType_changed(int geometryType){

    if (geometryType==0){ // simple geometry, flat bed
//       ui->length->setEnabled(true);ui->width->setEnabled(true);ui->depth->setEnabled(true);
       ui->customGridWidget->setVisible(false);

    } else { // Advanced geometry type

//        ui->length->setEnabled(false);ui->depth->setEnabled(false);ui->width->setEnabled(false);
        ui->customGridWidget->setVisible(true);

        int nGridPoints = ui->nGridPoints->value();
        QDoubleSpinBox *sp;


        if (nGridPoints>ui->geometry_table->columnCount()){
            for (int i =ui->geometry_table->columnCount();i<nGridPoints;i++){

                ui->geometry_table->insertColumn(i);

                for (int j = 0; j < 2 ; j++) {
                    ui->geometry_table->setCellWidget(j,i,new QDoubleSpinBox(ui->geometry_table));
                    sp = (QDoubleSpinBox*)ui->geometry_table->cellWidget(j,i);
                    sp->setDecimals(3);
                    sp->setMaximum(9999999);
                    if ((j==0) & (i == 0)){sp->setEnabled(false);}

                }
            }
        }

        if (nGridPoints<ui->geometry_table->columnCount()){
            for (int i=ui->geometry_table->columnCount();i>nGridPoints;i--){
                ui->geometry_table->removeColumn(i-1);
            }

        }


    }
}

// Call smooth functon CustomGrid::laplaceSmoothing
//
void MainWindow::smooth(){
        grid.laplaceSmoothing();

        // We always want a smooth grid, so we help the user a bit
        //
        for (int i=0;i<5;i++){
            grid.laplaceSmoothing();
        }

        // We update the grid file
        grid.writeGridToFile(ui->selectedGridFile->text());
}

// Read grid table and generate grid
//
void MainWindow::generateGrid(){
     grid.readTabel(ui->geometry_table,ui->nx->value(),ui->nz->value(),ui->sz->isChecked());
     grid.generateGrid();
     grid.writeGridToFile(ui->selectedGridFile->text());
     // Enable smooth and showGrid functions
     ui->showGrid->setEnabled(true);
     ui->smooth->setEnabled(true);
     ui->length->setValue(grid.L);
     ui->depth->setValue(grid.d);
}

void MainWindow::showGrid(){
    int plotType = ui->geometryType->currentIndex();
    if (plotType==0){
        grid.showGrid(ui->depth->value(),ui->length->value(),ui->nx->value(),ui->nz->value(),ui->sz->isChecked());
    } else if (plotType==1){
        grid.showGrid();
    }


}

void MainWindow::nGridPoints_changed(){

    geometryType_changed(1);
}

