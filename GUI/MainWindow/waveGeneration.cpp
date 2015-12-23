#include "mainwindow.h"
#include "ui_mainwindow.h"

void MainWindow::on_waveTheoryChanged()
{
    MainWindow::on_waveTheoryChanged(1);
}

void MainWindow::on_waveTheoryChanged(int waveTheory)
{
    if (waveTheory==1) { // stream function waves
        ui->widget_SF->setVisible(true);ui->LorP_ComboBox->setVisible(true);

        if (ui->LorP_ComboBox->currentIndex()==0){
            ui->SF_TL_unit->setText("m");
            ui->SF_TLabel->setText("Wave length");
        }
        if (ui->LorP_ComboBox->currentIndex()==1) {
            ui->SF_TL_unit->setText("s");
            ui->SF_TLabel->setText("Wave period");
        }
        ui->widget_JONSWAP->setVisible(false);
        ui->widget_waveFile->setVisible(false);
        ui->waveGeneration_widget->setEnabled(true);
        ui->widget_customSpectrum->setVisible(false);
    } else if (waveTheory==2){ // JONSWAP
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
        ui->widget_JONSWAP->setVisible(true);
        ui->widget_waveFile->setVisible(false);
        ui->waveGeneration_widget->setEnabled(true);
        ui->widget_customSpectrum->setVisible(false);
        ui->gamma_jonswap->setVisible(true);
        ui->gamma_jonswap_label->setVisible(true);
    } else if (waveTheory==3){
        ui->widget_JONSWAP->setVisible(true);
        ui->widget_SF->setVisible(false);
        ui->widget_waveFile->setVisible(false);
        ui->waveGeneration_widget->setEnabled(true);
        ui->widget_customSpectrum->setVisible(false);
        ui->gamma_jonswap->setVisible(false);
        ui->gamma_jonswap_label->setVisible(false);
    } else if (waveTheory==4) {
        ui->widget_customSpectrum->setVisible(true);
        ui->widget_waveFile->setVisible(false);
        ui->widget_JONSWAP->setVisible(false);
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
        ui->waveGeneration_widget->setEnabled(true);
    } else if (waveTheory==5) {
        ui->widget_waveFile->setVisible(true);
        ui->widget_JONSWAP->setVisible(false);
        ui->widget_customSpectrum->setVisible(false);
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
        ui->waveGeneration_widget->setEnabled(false);

    } else{
        ui->widget_customSpectrum->setVisible(false);
        ui->widget_JONSWAP->setVisible(false);
        ui->widget_SF->setVisible(false);ui->LorP_ComboBox->setVisible(false);
        ui->widget_waveFile->setVisible(false);
        ui->waveGeneration_widget->setEnabled(true);
    }
}

void MainWindow::selectWaveFile(){

    QString fileName = QFileInfo(QFileDialog::getOpenFileName(
                                     this,
                                     "Select file",
                                     dir.currentPath(),
                                     "All files (*.*);;"
                                     )).fileName();
    ui->selectedWaveFile->setText(fileName);


}

void MainWindow::selectWaveFile_eta(){

    QString fileName = QFileInfo(QFileDialog::getOpenFileName(
                                     this,
                                     "Select file",
                                     dir.currentPath(),
                                     "All files (*.*);;"
                                     )).fileName();
    ui->selectedWaveFile_eta->setText(fileName);


}
