/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created: Fri Apr 24 10:08:42 2015
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QPlainTextEdit>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTableWidget>
#include <QtGui/QTextBrowser>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionOpen;
    QWidget *centralWidget;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QTabWidget *output;
    QWidget *tab_11;
    QFrame *line_34;
    QLabel *label_79;
    QLabel *label_82;
    QFrame *line_35;
    QFrame *line_36;
    QLabel *label_87;
    QFrame *line_37;
    QPlainTextEdit *header_input;
    QLabel *label_89;
    QDoubleSpinBox *breaking_beta0;
    QLabel *label_90;
    QDoubleSpinBox *gravity_input;
    QLabel *label_91;
    QFrame *line_41;
    QLabel *label_92;
    QCheckBox *linear_onOff;
    QCheckBox *Nonlin_onOff;
    QFrame *line_65;
    QComboBox *initialCondition;
    QFrame *line_66;
    QLabel *label_139;
    QDoubleSpinBox *timeDuration;
    QLabel *label_141;
    QLabel *label_145;
    QDoubleSpinBox *tStart;
    QLabel *label_158;
    QDoubleSpinBox *Density;
    QLabel *label_159;
    QLabel *label_160;
    QLabel *label_161;
    QLabel *label_162;
    QLabel *label_165;
    QPushButton *OpenDirBrowser;
    QLineEdit *workingDir;
    QWidget *tab;
    QLabel *label_6;
    QFrame *line_3;
    QLabel *label_9;
    QLabel *label_5;
    QFrame *line_2;
    QFrame *line;
    QLabel *label_2;
    QCheckBox *sy;
    QFrame *line_4;
    QLabel *label_8;
    QLabel *label_10;
    QLabel *label_7;
    QCheckBox *sz;
    QCheckBox *sx;
    QLabel *label_4;
    QLabel *label_3;
    QLabel *label_88;
    QDoubleSpinBox *length;
    QDoubleSpinBox *width;
    QDoubleSpinBox *depth;
    QSpinBox *nx;
    QSpinBox *ny;
    QSpinBox *nz;
    QWidget *tab_2;
    QLabel *label_11;
    QLabel *label_12;
    QFrame *line_5;
    QFrame *line_6;
    QFrame *line_7;
    QFrame *line_8;
    QLabel *label_13;
    QLabel *label_27;
    QSpinBox *gamma;
    QSpinBox *alpha;
    QLabel *label_28;
    QLabel *label_29;
    QSpinBox *beta;
    QSpinBox *b;
    QLabel *label_30;
    QSpinBox *c;
    QLabel *label_31;
    QLabel *label_32;
    QSpinBox *a;
    QLabel *label_137;
    QDoubleSpinBox *dt;
    QLabel *label_138;
    QWidget *tab_5;
    QFrame *line_25;
    QLabel *label_52;
    QLabel *label_55;
    QFrame *line_27;
    QDoubleSpinBox *xGenStart;
    QDoubleSpinBox *xGenEnd;
    QDoubleSpinBox *yGenStart;
    QDoubleSpinBox *yGenEnd;
    QLabel *label_126;
    QLabel *label_127;
    QLabel *label_129;
    QLabel *label_130;
    QDoubleSpinBox *yAbsorbStart;
    QLabel *label_132;
    QDoubleSpinBox *xAbsorbStart;
    QLabel *label_133;
    QLabel *label_134;
    QDoubleSpinBox *xAbsorbEnd;
    QDoubleSpinBox *yAbsorbEnd;
    QLabel *label_135;
    QLabel *label_136;
    QFrame *line_26;
    QComboBox *waveType;
    QWidget *widget_SF;
    QLabel *SF_TLabel;
    QDoubleSpinBox *SF_H;
    QLabel *label_142;
    QDoubleSpinBox *SF_h;
    QDoubleSpinBox *SF_T;
    QLabel *label_140;
    QDoubleSpinBox *SF_U;
    QLabel *label_143;
    QLabel *label_144;
    QLabel *SF_TL_unit;
    QLabel *label_146;
    QLabel *label_147;
    QLabel *label_148;
    QComboBox *stokesOrEuler;
    QLabel *label_149;
    QSpinBox *SF_n;
    QLineEdit *irrFileName;
    QLabel *irrFileName_label;
    QComboBox *LorP_ComboBox;
    QComboBox *pressureDampingOrRelax;
    QWidget *widget_JONSWAP;
    QLabel *label_150;
    QDoubleSpinBox *Hs;
    QLabel *label_151;
    QDoubleSpinBox *h;
    QDoubleSpinBox *Tp;
    QLabel *label_152;
    QDoubleSpinBox *gamma_2;
    QLabel *label_153;
    QLabel *label_154;
    QLabel *label_155;
    QLabel *label_156;
    QLabel *label_157;
    QDoubleSpinBox *maxkh;
    QLabel *label_166;
    QLabel *label_167;
    QSpinBox *seed;
    QDoubleSpinBox *rampTime;
    QLabel *label_163;
    QLabel *label_164;
    QWidget *tab_9;
    QSpinBox *nOutFiles;
    QLabel *label;
    QFrame *line_28;
    QTableWidget *tableWidget;
    QCheckBox *storeAscii_onOff;
    QSpinBox *nASCII;
    QLabel *ACCII_label;
    QLabel *ASCII_label2;
    QWidget *tab_3;
    QPushButton *selectPPfile;
    QLineEdit *selectedPPfiles;
    QFrame *line_9;
    QPushButton *convert;
    QPushButton *selectGnuplotFile;
    QLineEdit *gnuplotFile;
    QLabel *label_14;
    QLabel *label_15;
    QPushButton *plot;
    QPushButton *read_bottom;
    QProgressBar *readProgressBar;
    QComboBox *SelectOutput;
    QWidget *morison_widget;
    QDoubleSpinBox *morison_x0;
    QLabel *label_16;
    QLabel *label_17;
    QDoubleSpinBox *morison_D;
    QLabel *label_18;
    QDoubleSpinBox *morison_cd;
    QLabel *label_19;
    QDoubleSpinBox *morison_cm;
    QLabel *convertStatus;
    QPlainTextEdit *convertWarning;
    QWidget *About;
    QComboBox *about_combobox;
    QTextBrowser *aboutText_OCW3D;
    QTextBrowser *aboutText_OCW3dGUI;
    QTextBrowser *aboutText_OCW3D_publications;
    QPushButton *run;
    QStatusBar *statusBar;
    QMenuBar *menuBar;
    QMenu *menuFile;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(816, 528);
        actionOpen = new QAction(MainWindow);
        actionOpen->setObjectName(QString::fromUtf8("actionOpen"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        pushButton = new QPushButton(centralWidget);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));
        pushButton->setGeometry(QRect(700, 430, 80, 23));
        pushButton_2 = new QPushButton(centralWidget);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));
        pushButton_2->setGeometry(QRect(549, 430, 131, 23));
        output = new QTabWidget(centralWidget);
        output->setObjectName(QString::fromUtf8("output"));
        output->setGeometry(QRect(0, 0, 771, 411));
        output->setLayoutDirection(Qt::LeftToRight);
        output->setTabPosition(QTabWidget::North);
        output->setTabShape(QTabWidget::Rounded);
        tab_11 = new QWidget();
        tab_11->setObjectName(QString::fromUtf8("tab_11"));
        line_34 = new QFrame(tab_11);
        line_34->setObjectName(QString::fromUtf8("line_34"));
        line_34->setGeometry(QRect(0, 160, 781, 16));
        line_34->setFrameShape(QFrame::HLine);
        line_34->setFrameShadow(QFrame::Sunken);
        label_79 = new QLabel(tab_11);
        label_79->setObjectName(QString::fromUtf8("label_79"));
        label_79->setGeometry(QRect(10, 10, 201, 16));
        label_82 = new QLabel(tab_11);
        label_82->setObjectName(QString::fromUtf8("label_82"));
        label_82->setGeometry(QRect(10, 70, 151, 16));
        line_35 = new QFrame(tab_11);
        line_35->setObjectName(QString::fromUtf8("line_35"));
        line_35->setGeometry(QRect(273, 70, 20, 91));
        line_35->setFrameShape(QFrame::VLine);
        line_35->setFrameShadow(QFrame::Sunken);
        line_36 = new QFrame(tab_11);
        line_36->setObjectName(QString::fromUtf8("line_36"));
        line_36->setGeometry(QRect(493, 70, 20, 91));
        line_36->setFrameShape(QFrame::VLine);
        line_36->setFrameShadow(QFrame::Sunken);
        label_87 = new QLabel(tab_11);
        label_87->setObjectName(QString::fromUtf8("label_87"));
        label_87->setGeometry(QRect(300, 70, 151, 16));
        line_37 = new QFrame(tab_11);
        line_37->setObjectName(QString::fromUtf8("line_37"));
        line_37->setGeometry(QRect(10, 60, 781, 16));
        line_37->setFrameShape(QFrame::HLine);
        line_37->setFrameShadow(QFrame::Sunken);
        header_input = new QPlainTextEdit(tab_11);
        header_input->setObjectName(QString::fromUtf8("header_input"));
        header_input->setGeometry(QRect(60, 10, 341, 31));
        label_89 = new QLabel(tab_11);
        label_89->setObjectName(QString::fromUtf8("label_89"));
        label_89->setGeometry(QRect(300, 100, 151, 16));
        breaking_beta0 = new QDoubleSpinBox(tab_11);
        breaking_beta0->setObjectName(QString::fromUtf8("breaking_beta0"));
        breaking_beta0->setGeometry(QRect(300, 120, 66, 24));
        breaking_beta0->setSingleStep(0.1);
        breaking_beta0->setValue(0);
        label_90 = new QLabel(tab_11);
        label_90->setObjectName(QString::fromUtf8("label_90"));
        label_90->setGeometry(QRect(520, 70, 151, 16));
        gravity_input = new QDoubleSpinBox(tab_11);
        gravity_input->setObjectName(QString::fromUtf8("gravity_input"));
        gravity_input->setGeometry(QRect(520, 120, 66, 24));
        gravity_input->setSingleStep(0.1);
        gravity_input->setValue(9.82);
        label_91 = new QLabel(tab_11);
        label_91->setObjectName(QString::fromUtf8("label_91"));
        label_91->setGeometry(QRect(520, 100, 51, 16));
        line_41 = new QFrame(tab_11);
        line_41->setObjectName(QString::fromUtf8("line_41"));
        line_41->setGeometry(QRect(493, 0, 20, 61));
        line_41->setFrameShape(QFrame::VLine);
        line_41->setFrameShadow(QFrame::Sunken);
        label_92 = new QLabel(tab_11);
        label_92->setObjectName(QString::fromUtf8("label_92"));
        label_92->setGeometry(QRect(520, 0, 151, 16));
        linear_onOff = new QCheckBox(tab_11);
        linear_onOff->setObjectName(QString::fromUtf8("linear_onOff"));
        linear_onOff->setGeometry(QRect(610, 20, 85, 21));
        Nonlin_onOff = new QCheckBox(tab_11);
        Nonlin_onOff->setObjectName(QString::fromUtf8("Nonlin_onOff"));
        Nonlin_onOff->setGeometry(QRect(520, 20, 85, 21));
        Nonlin_onOff->setCheckable(true);
        Nonlin_onOff->setChecked(true);
        Nonlin_onOff->setAutoExclusive(false);
        line_65 = new QFrame(tab_11);
        line_65->setObjectName(QString::fromUtf8("line_65"));
        line_65->setGeometry(QRect(273, 180, 20, 61));
        line_65->setFrameShape(QFrame::VLine);
        line_65->setFrameShadow(QFrame::Sunken);
        initialCondition = new QComboBox(tab_11);
        initialCondition->setObjectName(QString::fromUtf8("initialCondition"));
        initialCondition->setGeometry(QRect(10, 110, 121, 23));
        line_66 = new QFrame(tab_11);
        line_66->setObjectName(QString::fromUtf8("line_66"));
        line_66->setGeometry(QRect(0, 240, 781, 16));
        line_66->setFrameShape(QFrame::HLine);
        line_66->setFrameShadow(QFrame::Sunken);
        label_139 = new QLabel(tab_11);
        label_139->setObjectName(QString::fromUtf8("label_139"));
        label_139->setGeometry(QRect(10, 180, 31, 16));
        timeDuration = new QDoubleSpinBox(tab_11);
        timeDuration->setObjectName(QString::fromUtf8("timeDuration"));
        timeDuration->setGeometry(QRect(10, 220, 71, 24));
        timeDuration->setDecimals(2);
        timeDuration->setMaximum(100000);
        timeDuration->setSingleStep(0.1);
        timeDuration->setValue(10800);
        label_141 = new QLabel(tab_11);
        label_141->setObjectName(QString::fromUtf8("label_141"));
        label_141->setGeometry(QRect(10, 200, 71, 16));
        label_145 = new QLabel(tab_11);
        label_145->setObjectName(QString::fromUtf8("label_145"));
        label_145->setGeometry(QRect(110, 200, 71, 16));
        tStart = new QDoubleSpinBox(tab_11);
        tStart->setObjectName(QString::fromUtf8("tStart"));
        tStart->setGeometry(QRect(110, 220, 71, 24));
        tStart->setDecimals(2);
        tStart->setMaximum(100000);
        tStart->setSingleStep(0.1);
        tStart->setValue(0);
        label_158 = new QLabel(tab_11);
        label_158->setObjectName(QString::fromUtf8("label_158"));
        label_158->setGeometry(QRect(640, 100, 51, 16));
        Density = new QDoubleSpinBox(tab_11);
        Density->setObjectName(QString::fromUtf8("Density"));
        Density->setGeometry(QRect(640, 120, 66, 24));
        Density->setMaximum(5000);
        Density->setSingleStep(0.1);
        Density->setValue(1025);
        label_159 = new QLabel(tab_11);
        label_159->setObjectName(QString::fromUtf8("label_159"));
        label_159->setGeometry(QRect(590, 120, 51, 16));
        label_160 = new QLabel(tab_11);
        label_160->setObjectName(QString::fromUtf8("label_160"));
        label_160->setGeometry(QRect(710, 120, 51, 16));
        label_161 = new QLabel(tab_11);
        label_161->setObjectName(QString::fromUtf8("label_161"));
        label_161->setGeometry(QRect(190, 220, 16, 16));
        label_162 = new QLabel(tab_11);
        label_162->setObjectName(QString::fromUtf8("label_162"));
        label_162->setGeometry(QRect(90, 220, 16, 16));
        label_165 = new QLabel(tab_11);
        label_165->setObjectName(QString::fromUtf8("label_165"));
        label_165->setGeometry(QRect(300, 180, 111, 16));
        OpenDirBrowser = new QPushButton(tab_11);
        OpenDirBrowser->setObjectName(QString::fromUtf8("OpenDirBrowser"));
        OpenDirBrowser->setGeometry(QRect(570, 200, 91, 24));
        workingDir = new QLineEdit(tab_11);
        workingDir->setObjectName(QString::fromUtf8("workingDir"));
        workingDir->setEnabled(false);
        workingDir->setGeometry(QRect(300, 200, 241, 23));
        output->addTab(tab_11, QString());
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        label_6 = new QLabel(tab);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(290, 80, 57, 15));
        line_3 = new QFrame(tab);
        line_3->setObjectName(QString::fromUtf8("line_3"));
        line_3->setGeometry(QRect(-10, 120, 781, 16));
        line_3->setFrameShape(QFrame::HLine);
        line_3->setFrameShadow(QFrame::Sunken);
        label_9 = new QLabel(tab);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(550, 60, 71, 16));
        label_5 = new QLabel(tab);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(170, 80, 75, 16));
        label_5->setMinimumSize(QSize(75, 0));
        line_2 = new QFrame(tab);
        line_2->setObjectName(QString::fromUtf8("line_2"));
        line_2->setGeometry(QRect(280, 60, 3, 61));
        line_2->setFrameShape(QFrame::VLine);
        line_2->setFrameShadow(QFrame::Sunken);
        line = new QFrame(tab);
        line->setObjectName(QString::fromUtf8("line"));
        line->setGeometry(QRect(0, 50, 781, 16));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);
        label_2 = new QLabel(tab);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(10, 80, 75, 16));
        label_2->setMinimumSize(QSize(75, 0));
        sy = new QCheckBox(tab);
        sy->setObjectName(QString::fromUtf8("sy"));
        sy->setGeometry(QRect(590, 100, 31, 21));
        sy->setChecked(false);
        line_4 = new QFrame(tab);
        line_4->setObjectName(QString::fromUtf8("line_4"));
        line_4->setGeometry(QRect(540, 60, 3, 61));
        line_4->setFrameShape(QFrame::VLine);
        line_4->setFrameShadow(QFrame::Sunken);
        label_8 = new QLabel(tab);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(450, 80, 57, 15));
        label_10 = new QLabel(tab);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setGeometry(QRect(290, 60, 71, 16));
        label_7 = new QLabel(tab);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(370, 80, 57, 15));
        sz = new QCheckBox(tab);
        sz->setObjectName(QString::fromUtf8("sz"));
        sz->setGeometry(QRect(630, 100, 41, 21));
        sz->setChecked(true);
        sx = new QCheckBox(tab);
        sx->setObjectName(QString::fromUtf8("sx"));
        sx->setGeometry(QRect(550, 100, 31, 21));
        sx->setChecked(false);
        label_4 = new QLabel(tab);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(10, 60, 57, 15));
        label_3 = new QLabel(tab);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(90, 80, 75, 16));
        label_3->setMinimumSize(QSize(75, 0));
        label_88 = new QLabel(tab);
        label_88->setObjectName(QString::fromUtf8("label_88"));
        label_88->setGeometry(QRect(10, 10, 201, 16));
        length = new QDoubleSpinBox(tab);
        length->setObjectName(QString::fromUtf8("length"));
        length->setGeometry(QRect(10, 100, 66, 24));
        length->setMaximum(10000);
        length->setSingleStep(0.1);
        width = new QDoubleSpinBox(tab);
        width->setObjectName(QString::fromUtf8("width"));
        width->setGeometry(QRect(90, 100, 66, 24));
        width->setMaximum(10000);
        width->setSingleStep(0.1);
        depth = new QDoubleSpinBox(tab);
        depth->setObjectName(QString::fromUtf8("depth"));
        depth->setGeometry(QRect(170, 100, 66, 24));
        depth->setMaximum(10000);
        depth->setSingleStep(0.1);
        nx = new QSpinBox(tab);
        nx->setObjectName(QString::fromUtf8("nx"));
        nx->setGeometry(QRect(290, 100, 61, 24));
        nx->setMaximum(10000);
        ny = new QSpinBox(tab);
        ny->setObjectName(QString::fromUtf8("ny"));
        ny->setGeometry(QRect(370, 100, 61, 24));
        ny->setMaximum(10000);
        nz = new QSpinBox(tab);
        nz->setObjectName(QString::fromUtf8("nz"));
        nz->setGeometry(QRect(450, 100, 61, 24));
        nz->setMaximum(10000);
        output->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        label_11 = new QLabel(tab_2);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setGeometry(QRect(10, 10, 201, 16));
        label_12 = new QLabel(tab_2);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setGeometry(QRect(300, 60, 151, 16));
        line_5 = new QFrame(tab_2);
        line_5->setObjectName(QString::fromUtf8("line_5"));
        line_5->setGeometry(QRect(0, 120, 781, 16));
        line_5->setFrameShape(QFrame::HLine);
        line_5->setFrameShadow(QFrame::Sunken);
        line_6 = new QFrame(tab_2);
        line_6->setObjectName(QString::fromUtf8("line_6"));
        line_6->setGeometry(QRect(290, 60, 3, 61));
        line_6->setFrameShape(QFrame::VLine);
        line_6->setFrameShadow(QFrame::Sunken);
        line_7 = new QFrame(tab_2);
        line_7->setObjectName(QString::fromUtf8("line_7"));
        line_7->setGeometry(QRect(10, 50, 781, 16));
        line_7->setFrameShape(QFrame::HLine);
        line_7->setFrameShadow(QFrame::Sunken);
        line_8 = new QFrame(tab_2);
        line_8->setObjectName(QString::fromUtf8("line_8"));
        line_8->setGeometry(QRect(510, 60, 3, 61));
        line_8->setFrameShape(QFrame::VLine);
        line_8->setFrameShadow(QFrame::Sunken);
        label_13 = new QLabel(tab_2);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        label_13->setGeometry(QRect(10, 60, 151, 16));
        label_27 = new QLabel(tab_2);
        label_27->setObjectName(QString::fromUtf8("label_27"));
        label_27->setGeometry(QRect(360, 80, 57, 15));
        gamma = new QSpinBox(tab_2);
        gamma->setObjectName(QString::fromUtf8("gamma"));
        gamma->setGeometry(QRect(169, 100, 61, 24));
        gamma->setMinimumSize(QSize(50, 0));
        gamma->setMaximum(10);
        gamma->setValue(3);
        alpha = new QSpinBox(tab_2);
        alpha->setObjectName(QString::fromUtf8("alpha"));
        alpha->setGeometry(QRect(9, 100, 61, 24));
        alpha->setMinimumSize(QSize(50, 0));
        alpha->setMaximum(10);
        alpha->setValue(3);
        label_28 = new QLabel(tab_2);
        label_28->setObjectName(QString::fromUtf8("label_28"));
        label_28->setGeometry(QRect(170, 80, 75, 16));
        label_28->setMinimumSize(QSize(75, 0));
        label_29 = new QLabel(tab_2);
        label_29->setObjectName(QString::fromUtf8("label_29"));
        label_29->setGeometry(QRect(90, 80, 75, 16));
        label_29->setMinimumSize(QSize(75, 0));
        beta = new QSpinBox(tab_2);
        beta->setObjectName(QString::fromUtf8("beta"));
        beta->setGeometry(QRect(90, 100, 61, 24));
        beta->setMaximum(10);
        beta->setValue(3);
        b = new QSpinBox(tab_2);
        b->setObjectName(QString::fromUtf8("b"));
        b->setGeometry(QRect(360, 100, 51, 24));
        b->setMaximum(10);
        b->setValue(1);
        label_30 = new QLabel(tab_2);
        label_30->setObjectName(QString::fromUtf8("label_30"));
        label_30->setGeometry(QRect(300, 80, 57, 15));
        c = new QSpinBox(tab_2);
        c->setObjectName(QString::fromUtf8("c"));
        c->setGeometry(QRect(420, 100, 51, 24));
        c->setMaximum(10);
        c->setValue(1);
        label_31 = new QLabel(tab_2);
        label_31->setObjectName(QString::fromUtf8("label_31"));
        label_31->setGeometry(QRect(420, 80, 57, 15));
        label_32 = new QLabel(tab_2);
        label_32->setObjectName(QString::fromUtf8("label_32"));
        label_32->setGeometry(QRect(10, 80, 75, 16));
        label_32->setMinimumSize(QSize(75, 0));
        a = new QSpinBox(tab_2);
        a->setObjectName(QString::fromUtf8("a"));
        a->setGeometry(QRect(300, 100, 51, 24));
        a->setMaximum(10);
        a->setValue(1);
        label_137 = new QLabel(tab_2);
        label_137->setObjectName(QString::fromUtf8("label_137"));
        label_137->setGeometry(QRect(520, 60, 151, 16));
        dt = new QDoubleSpinBox(tab_2);
        dt->setObjectName(QString::fromUtf8("dt"));
        dt->setGeometry(QRect(520, 100, 66, 24));
        dt->setValue(0.01);
        label_138 = new QLabel(tab_2);
        label_138->setObjectName(QString::fromUtf8("label_138"));
        label_138->setGeometry(QRect(520, 80, 57, 15));
        output->addTab(tab_2, QString());
        tab_5 = new QWidget();
        tab_5->setObjectName(QString::fromUtf8("tab_5"));
        line_25 = new QFrame(tab_5);
        line_25->setObjectName(QString::fromUtf8("line_25"));
        line_25->setGeometry(QRect(0, 270, 781, 16));
        line_25->setFrameShape(QFrame::HLine);
        line_25->setFrameShadow(QFrame::Sunken);
        label_52 = new QLabel(tab_5);
        label_52->setObjectName(QString::fromUtf8("label_52"));
        label_52->setGeometry(QRect(10, 20, 331, 16));
        label_55 = new QLabel(tab_5);
        label_55->setObjectName(QString::fromUtf8("label_55"));
        label_55->setGeometry(QRect(10, 140, 151, 16));
        line_27 = new QFrame(tab_5);
        line_27->setObjectName(QString::fromUtf8("line_27"));
        line_27->setGeometry(QRect(10, 130, 781, 16));
        line_27->setFrameShape(QFrame::HLine);
        line_27->setFrameShadow(QFrame::Sunken);
        xGenStart = new QDoubleSpinBox(tab_5);
        xGenStart->setObjectName(QString::fromUtf8("xGenStart"));
        xGenStart->setGeometry(QRect(50, 190, 66, 24));
        xGenStart->setMaximum(100000);
        xGenEnd = new QDoubleSpinBox(tab_5);
        xGenEnd->setObjectName(QString::fromUtf8("xGenEnd"));
        xGenEnd->setGeometry(QRect(50, 220, 66, 24));
        xGenEnd->setMaximum(100000);
        yGenStart = new QDoubleSpinBox(tab_5);
        yGenStart->setObjectName(QString::fromUtf8("yGenStart"));
        yGenStart->setGeometry(QRect(140, 190, 66, 24));
        yGenStart->setMaximum(100000);
        yGenEnd = new QDoubleSpinBox(tab_5);
        yGenEnd->setObjectName(QString::fromUtf8("yGenEnd"));
        yGenEnd->setGeometry(QRect(140, 220, 66, 24));
        yGenEnd->setMaximum(100000);
        label_126 = new QLabel(tab_5);
        label_126->setObjectName(QString::fromUtf8("label_126"));
        label_126->setGeometry(QRect(50, 170, 41, 16));
        label_127 = new QLabel(tab_5);
        label_127->setObjectName(QString::fromUtf8("label_127"));
        label_127->setGeometry(QRect(140, 170, 41, 16));
        label_129 = new QLabel(tab_5);
        label_129->setObjectName(QString::fromUtf8("label_129"));
        label_129->setGeometry(QRect(10, 190, 41, 16));
        label_130 = new QLabel(tab_5);
        label_130->setObjectName(QString::fromUtf8("label_130"));
        label_130->setGeometry(QRect(10, 220, 31, 16));
        yAbsorbStart = new QDoubleSpinBox(tab_5);
        yAbsorbStart->setObjectName(QString::fromUtf8("yAbsorbStart"));
        yAbsorbStart->setGeometry(QRect(530, 190, 66, 24));
        yAbsorbStart->setMaximum(100000);
        label_132 = new QLabel(tab_5);
        label_132->setObjectName(QString::fromUtf8("label_132"));
        label_132->setGeometry(QRect(440, 170, 41, 16));
        xAbsorbStart = new QDoubleSpinBox(tab_5);
        xAbsorbStart->setObjectName(QString::fromUtf8("xAbsorbStart"));
        xAbsorbStart->setGeometry(QRect(440, 190, 66, 24));
        xAbsorbStart->setMaximum(100000);
        label_133 = new QLabel(tab_5);
        label_133->setObjectName(QString::fromUtf8("label_133"));
        label_133->setGeometry(QRect(400, 220, 31, 16));
        label_134 = new QLabel(tab_5);
        label_134->setObjectName(QString::fromUtf8("label_134"));
        label_134->setGeometry(QRect(400, 140, 151, 16));
        xAbsorbEnd = new QDoubleSpinBox(tab_5);
        xAbsorbEnd->setObjectName(QString::fromUtf8("xAbsorbEnd"));
        xAbsorbEnd->setGeometry(QRect(440, 220, 66, 24));
        xAbsorbEnd->setMaximum(100000);
        yAbsorbEnd = new QDoubleSpinBox(tab_5);
        yAbsorbEnd->setObjectName(QString::fromUtf8("yAbsorbEnd"));
        yAbsorbEnd->setGeometry(QRect(530, 220, 66, 24));
        yAbsorbEnd->setMaximum(100000);
        label_135 = new QLabel(tab_5);
        label_135->setObjectName(QString::fromUtf8("label_135"));
        label_135->setGeometry(QRect(530, 170, 41, 16));
        label_136 = new QLabel(tab_5);
        label_136->setObjectName(QString::fromUtf8("label_136"));
        label_136->setGeometry(QRect(400, 190, 41, 16));
        line_26 = new QFrame(tab_5);
        line_26->setObjectName(QString::fromUtf8("line_26"));
        line_26->setGeometry(QRect(363, 150, 20, 121));
        line_26->setFrameShape(QFrame::VLine);
        line_26->setFrameShadow(QFrame::Sunken);
        waveType = new QComboBox(tab_5);
        waveType->setObjectName(QString::fromUtf8("waveType"));
        waveType->setGeometry(QRect(250, 20, 121, 23));
        widget_SF = new QWidget(tab_5);
        widget_SF->setObjectName(QString::fromUtf8("widget_SF"));
        widget_SF->setEnabled(true);
        widget_SF->setGeometry(QRect(-40, 60, 771, 80));
        SF_TLabel = new QLabel(widget_SF);
        SF_TLabel->setObjectName(QString::fromUtf8("SF_TLabel"));
        SF_TLabel->setGeometry(QRect(150, 10, 81, 16));
        SF_H = new QDoubleSpinBox(widget_SF);
        SF_H->setObjectName(QString::fromUtf8("SF_H"));
        SF_H->setGeometry(QRect(50, 30, 66, 24));
        label_142 = new QLabel(widget_SF);
        label_142->setObjectName(QString::fromUtf8("label_142"));
        label_142->setGeometry(QRect(250, 10, 81, 16));
        SF_h = new QDoubleSpinBox(widget_SF);
        SF_h->setObjectName(QString::fromUtf8("SF_h"));
        SF_h->setGeometry(QRect(250, 30, 66, 24));
        SF_h->setMaximum(9999);
        SF_T = new QDoubleSpinBox(widget_SF);
        SF_T->setObjectName(QString::fromUtf8("SF_T"));
        SF_T->setGeometry(QRect(150, 30, 66, 24));
        SF_T->setMaximum(1000);
        SF_T->setSingleStep(0.1);
        label_140 = new QLabel(widget_SF);
        label_140->setObjectName(QString::fromUtf8("label_140"));
        label_140->setGeometry(QRect(50, 10, 81, 16));
        SF_U = new QDoubleSpinBox(widget_SF);
        SF_U->setObjectName(QString::fromUtf8("SF_U"));
        SF_U->setGeometry(QRect(350, 30, 66, 24));
        label_143 = new QLabel(widget_SF);
        label_143->setObjectName(QString::fromUtf8("label_143"));
        label_143->setGeometry(QRect(350, 10, 101, 16));
        label_144 = new QLabel(widget_SF);
        label_144->setObjectName(QString::fromUtf8("label_144"));
        label_144->setGeometry(QRect(120, 30, 21, 16));
        SF_TL_unit = new QLabel(widget_SF);
        SF_TL_unit->setObjectName(QString::fromUtf8("SF_TL_unit"));
        SF_TL_unit->setGeometry(QRect(220, 30, 21, 16));
        label_146 = new QLabel(widget_SF);
        label_146->setObjectName(QString::fromUtf8("label_146"));
        label_146->setGeometry(QRect(320, 30, 21, 16));
        label_147 = new QLabel(widget_SF);
        label_147->setObjectName(QString::fromUtf8("label_147"));
        label_147->setGeometry(QRect(420, 30, 21, 16));
        label_148 = new QLabel(widget_SF);
        label_148->setObjectName(QString::fromUtf8("label_148"));
        label_148->setGeometry(QRect(460, 10, 81, 16));
        stokesOrEuler = new QComboBox(widget_SF);
        stokesOrEuler->setObjectName(QString::fromUtf8("stokesOrEuler"));
        stokesOrEuler->setGeometry(QRect(460, 30, 91, 23));
        label_149 = new QLabel(widget_SF);
        label_149->setObjectName(QString::fromUtf8("label_149"));
        label_149->setGeometry(QRect(570, 10, 81, 16));
        SF_n = new QSpinBox(widget_SF);
        SF_n->setObjectName(QString::fromUtf8("SF_n"));
        SF_n->setGeometry(QRect(570, 30, 47, 24));
        SF_n->setMinimum(1);
        SF_n->setMaximum(40);
        SF_n->setValue(32);
        irrFileName = new QLineEdit(widget_SF);
        irrFileName->setObjectName(QString::fromUtf8("irrFileName"));
        irrFileName->setEnabled(false);
        irrFileName->setGeometry(QRect(640, 30, 113, 23));
        irrFileName_label = new QLabel(widget_SF);
        irrFileName_label->setObjectName(QString::fromUtf8("irrFileName_label"));
        irrFileName_label->setGeometry(QRect(640, 10, 81, 16));
        SF_TLabel->raise();
        SF_T->raise();
        SF_h->raise();
        label_140->raise();
        SF_H->raise();
        label_142->raise();
        SF_U->raise();
        label_143->raise();
        label_144->raise();
        SF_TL_unit->raise();
        label_146->raise();
        label_147->raise();
        label_148->raise();
        stokesOrEuler->raise();
        label_149->raise();
        SF_n->raise();
        irrFileName->raise();
        irrFileName_label->raise();
        LorP_ComboBox = new QComboBox(tab_5);
        LorP_ComboBox->setObjectName(QString::fromUtf8("LorP_ComboBox"));
        LorP_ComboBox->setGeometry(QRect(410, 20, 111, 23));
        pressureDampingOrRelax = new QComboBox(tab_5);
        pressureDampingOrRelax->setObjectName(QString::fromUtf8("pressureDampingOrRelax"));
        pressureDampingOrRelax->setGeometry(QRect(610, 190, 141, 23));
        widget_JONSWAP = new QWidget(tab_5);
        widget_JONSWAP->setObjectName(QString::fromUtf8("widget_JONSWAP"));
        widget_JONSWAP->setEnabled(true);
        widget_JONSWAP->setGeometry(QRect(-10, 290, 681, 80));
        label_150 = new QLabel(widget_JONSWAP);
        label_150->setObjectName(QString::fromUtf8("label_150"));
        label_150->setGeometry(QRect(150, 10, 81, 16));
        Hs = new QDoubleSpinBox(widget_JONSWAP);
        Hs->setObjectName(QString::fromUtf8("Hs"));
        Hs->setGeometry(QRect(50, 30, 66, 24));
        label_151 = new QLabel(widget_JONSWAP);
        label_151->setObjectName(QString::fromUtf8("label_151"));
        label_151->setGeometry(QRect(250, 10, 81, 16));
        h = new QDoubleSpinBox(widget_JONSWAP);
        h->setObjectName(QString::fromUtf8("h"));
        h->setGeometry(QRect(250, 30, 66, 24));
        Tp = new QDoubleSpinBox(widget_JONSWAP);
        Tp->setObjectName(QString::fromUtf8("Tp"));
        Tp->setGeometry(QRect(150, 30, 66, 24));
        label_152 = new QLabel(widget_JONSWAP);
        label_152->setObjectName(QString::fromUtf8("label_152"));
        label_152->setGeometry(QRect(50, 10, 81, 16));
        gamma_2 = new QDoubleSpinBox(widget_JONSWAP);
        gamma_2->setObjectName(QString::fromUtf8("gamma_2"));
        gamma_2->setGeometry(QRect(350, 30, 66, 24));
        gamma_2->setValue(3.3);
        label_153 = new QLabel(widget_JONSWAP);
        label_153->setObjectName(QString::fromUtf8("label_153"));
        label_153->setGeometry(QRect(350, 10, 81, 16));
        label_154 = new QLabel(widget_JONSWAP);
        label_154->setObjectName(QString::fromUtf8("label_154"));
        label_154->setGeometry(QRect(120, 30, 21, 16));
        label_155 = new QLabel(widget_JONSWAP);
        label_155->setObjectName(QString::fromUtf8("label_155"));
        label_155->setGeometry(QRect(220, 30, 21, 16));
        label_156 = new QLabel(widget_JONSWAP);
        label_156->setObjectName(QString::fromUtf8("label_156"));
        label_156->setGeometry(QRect(320, 30, 21, 16));
        label_157 = new QLabel(widget_JONSWAP);
        label_157->setObjectName(QString::fromUtf8("label_157"));
        label_157->setGeometry(QRect(420, 30, 21, 16));
        maxkh = new QDoubleSpinBox(widget_JONSWAP);
        maxkh->setObjectName(QString::fromUtf8("maxkh"));
        maxkh->setGeometry(QRect(450, 30, 66, 24));
        maxkh->setValue(5);
        label_166 = new QLabel(widget_JONSWAP);
        label_166->setObjectName(QString::fromUtf8("label_166"));
        label_166->setGeometry(QRect(450, 10, 81, 16));
        label_167 = new QLabel(widget_JONSWAP);
        label_167->setObjectName(QString::fromUtf8("label_167"));
        label_167->setGeometry(QRect(540, 10, 81, 16));
        seed = new QSpinBox(widget_JONSWAP);
        seed->setObjectName(QString::fromUtf8("seed"));
        seed->setGeometry(QRect(540, 30, 47, 24));
        seed->setMinimum(1);
        seed->setMaximum(5);
        rampTime = new QDoubleSpinBox(tab_5);
        rampTime->setObjectName(QString::fromUtf8("rampTime"));
        rampTime->setGeometry(QRect(230, 190, 66, 24));
        label_163 = new QLabel(tab_5);
        label_163->setObjectName(QString::fromUtf8("label_163"));
        label_163->setGeometry(QRect(230, 170, 91, 16));
        label_164 = new QLabel(tab_5);
        label_164->setObjectName(QString::fromUtf8("label_164"));
        label_164->setGeometry(QRect(620, 170, 41, 16));
        output->addTab(tab_5, QString());
        widget_SF->raise();
        line_25->raise();
        label_52->raise();
        label_55->raise();
        line_27->raise();
        xGenStart->raise();
        xGenEnd->raise();
        yGenStart->raise();
        yGenEnd->raise();
        label_126->raise();
        label_127->raise();
        label_129->raise();
        label_130->raise();
        yAbsorbStart->raise();
        label_132->raise();
        xAbsorbStart->raise();
        label_133->raise();
        label_134->raise();
        xAbsorbEnd->raise();
        yAbsorbEnd->raise();
        label_135->raise();
        label_136->raise();
        line_26->raise();
        waveType->raise();
        LorP_ComboBox->raise();
        pressureDampingOrRelax->raise();
        widget_JONSWAP->raise();
        rampTime->raise();
        label_163->raise();
        label_164->raise();
        tab_9 = new QWidget();
        tab_9->setObjectName(QString::fromUtf8("tab_9"));
        nOutFiles = new QSpinBox(tab_9);
        nOutFiles->setObjectName(QString::fromUtf8("nOutFiles"));
        nOutFiles->setGeometry(QRect(160, 20, 47, 24));
        nOutFiles->setMaximum(9);
        label = new QLabel(tab_9);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(10, 20, 151, 16));
        line_28 = new QFrame(tab_9);
        line_28->setObjectName(QString::fromUtf8("line_28"));
        line_28->setGeometry(QRect(7, 50, 761, 20));
        line_28->setFrameShape(QFrame::HLine);
        line_28->setFrameShadow(QFrame::Sunken);
        tableWidget = new QTableWidget(tab_9);
        tableWidget->setObjectName(QString::fromUtf8("tableWidget"));
        tableWidget->setGeometry(QRect(45, 100, 621, 241));
        storeAscii_onOff = new QCheckBox(tab_9);
        storeAscii_onOff->setObjectName(QString::fromUtf8("storeAscii_onOff"));
        storeAscii_onOff->setGeometry(QRect(240, 20, 121, 21));
        nASCII = new QSpinBox(tab_9);
        nASCII->setObjectName(QString::fromUtf8("nASCII"));
        nASCII->setGeometry(QRect(570, 20, 61, 24));
        nASCII->setMaximum(1000);
        nASCII->setValue(20);
        ACCII_label = new QLabel(tab_9);
        ACCII_label->setObjectName(QString::fromUtf8("ACCII_label"));
        ACCII_label->setGeometry(QRect(380, 20, 211, 16));
        ASCII_label2 = new QLabel(tab_9);
        ASCII_label2->setObjectName(QString::fromUtf8("ASCII_label2"));
        ASCII_label2->setGeometry(QRect(640, 20, 211, 16));
        output->addTab(tab_9, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QString::fromUtf8("tab_3"));
        selectPPfile = new QPushButton(tab_3);
        selectPPfile->setObjectName(QString::fromUtf8("selectPPfile"));
        selectPPfile->setGeometry(QRect(10, 140, 91, 24));
        selectedPPfiles = new QLineEdit(tab_3);
        selectedPPfiles->setObjectName(QString::fromUtf8("selectedPPfiles"));
        selectedPPfiles->setGeometry(QRect(110, 140, 251, 23));
        line_9 = new QFrame(tab_3);
        line_9->setObjectName(QString::fromUtf8("line_9"));
        line_9->setGeometry(QRect(0, 80, 741, 16));
        line_9->setFrameShape(QFrame::HLine);
        line_9->setFrameShadow(QFrame::Sunken);
        convert = new QPushButton(tab_3);
        convert->setObjectName(QString::fromUtf8("convert"));
        convert->setEnabled(true);
        convert->setGeometry(QRect(140, 190, 91, 24));
        selectGnuplotFile = new QPushButton(tab_3);
        selectGnuplotFile->setObjectName(QString::fromUtf8("selectGnuplotFile"));
        selectGnuplotFile->setGeometry(QRect(10, 40, 91, 24));
        gnuplotFile = new QLineEdit(tab_3);
        gnuplotFile->setObjectName(QString::fromUtf8("gnuplotFile"));
        gnuplotFile->setGeometry(QRect(110, 40, 251, 23));
        label_14 = new QLabel(tab_3);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setGeometry(QRect(13, 10, 81, 21));
        label_15 = new QLabel(tab_3);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setGeometry(QRect(20, 100, 81, 21));
        plot = new QPushButton(tab_3);
        plot->setObjectName(QString::fromUtf8("plot"));
        plot->setGeometry(QRect(380, 40, 91, 24));
        read_bottom = new QPushButton(tab_3);
        read_bottom->setObjectName(QString::fromUtf8("read_bottom"));
        read_bottom->setGeometry(QRect(380, 140, 91, 24));
        readProgressBar = new QProgressBar(tab_3);
        readProgressBar->setObjectName(QString::fromUtf8("readProgressBar"));
        readProgressBar->setGeometry(QRect(490, 140, 141, 23));
        readProgressBar->setValue(24);
        readProgressBar->setTextVisible(true);
        SelectOutput = new QComboBox(tab_3);
        SelectOutput->setObjectName(QString::fromUtf8("SelectOutput"));
        SelectOutput->setGeometry(QRect(10, 190, 111, 24));
        morison_widget = new QWidget(tab_3);
        morison_widget->setObjectName(QString::fromUtf8("morison_widget"));
        morison_widget->setEnabled(true);
        morison_widget->setGeometry(QRect(20, 230, 351, 81));
        morison_x0 = new QDoubleSpinBox(morison_widget);
        morison_x0->setObjectName(QString::fromUtf8("morison_x0"));
        morison_x0->setGeometry(QRect(0, 40, 62, 23));
        morison_x0->setMaximum(1000);
        morison_x0->setValue(501);
        label_16 = new QLabel(morison_widget);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setGeometry(QRect(10, 10, 41, 21));
        label_17 = new QLabel(morison_widget);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        label_17->setGeometry(QRect(100, 10, 41, 21));
        morison_D = new QDoubleSpinBox(morison_widget);
        morison_D->setObjectName(QString::fromUtf8("morison_D"));
        morison_D->setGeometry(QRect(90, 40, 62, 23));
        morison_D->setSingleStep(0.1);
        morison_D->setValue(5);
        label_18 = new QLabel(morison_widget);
        label_18->setObjectName(QString::fromUtf8("label_18"));
        label_18->setGeometry(QRect(190, 10, 41, 21));
        morison_cd = new QDoubleSpinBox(morison_widget);
        morison_cd->setObjectName(QString::fromUtf8("morison_cd"));
        morison_cd->setGeometry(QRect(180, 40, 62, 23));
        morison_cd->setSingleStep(0.1);
        morison_cd->setValue(0.95);
        label_19 = new QLabel(morison_widget);
        label_19->setObjectName(QString::fromUtf8("label_19"));
        label_19->setGeometry(QRect(280, 10, 41, 21));
        morison_cm = new QDoubleSpinBox(morison_widget);
        morison_cm->setObjectName(QString::fromUtf8("morison_cm"));
        morison_cm->setGeometry(QRect(270, 40, 62, 23));
        morison_cm->setValue(2);
        convertStatus = new QLabel(tab_3);
        convertStatus->setObjectName(QString::fromUtf8("convertStatus"));
        convertStatus->setGeometry(QRect(260, 190, 141, 16));
        convertWarning = new QPlainTextEdit(tab_3);
        convertWarning->setObjectName(QString::fromUtf8("convertWarning"));
        convertWarning->setGeometry(QRect(480, 130, 231, 51));
        convertWarning->setAutoFillBackground(true);
        convertWarning->setReadOnly(true);
        convertWarning->setBackgroundVisible(false);
        output->addTab(tab_3, QString());
        About = new QWidget();
        About->setObjectName(QString::fromUtf8("About"));
        About->setLayoutDirection(Qt::RightToLeft);
        about_combobox = new QComboBox(About);
        about_combobox->setObjectName(QString::fromUtf8("about_combobox"));
        about_combobox->setGeometry(QRect(580, 20, 151, 24));
        about_combobox->setLayoutDirection(Qt::LeftToRight);
        aboutText_OCW3D = new QTextBrowser(About);
        aboutText_OCW3D->setObjectName(QString::fromUtf8("aboutText_OCW3D"));
        aboutText_OCW3D->setGeometry(QRect(0, 10, 541, 341));
        aboutText_OCW3D->setLayoutDirection(Qt::LeftToRight);
        aboutText_OCW3dGUI = new QTextBrowser(About);
        aboutText_OCW3dGUI->setObjectName(QString::fromUtf8("aboutText_OCW3dGUI"));
        aboutText_OCW3dGUI->setGeometry(QRect(20, 50, 411, 111));
        aboutText_OCW3dGUI->setLayoutDirection(Qt::LeftToRight);
        aboutText_OCW3D_publications = new QTextBrowser(About);
        aboutText_OCW3D_publications->setObjectName(QString::fromUtf8("aboutText_OCW3D_publications"));
        aboutText_OCW3D_publications->setGeometry(QRect(80, 20, 441, 461));
        aboutText_OCW3D_publications->setLayoutDirection(Qt::LeftToRight);
        output->addTab(About, QString());
        run = new QPushButton(centralWidget);
        run->setObjectName(QString::fromUtf8("run"));
        run->setEnabled(false);
        run->setGeometry(QRect(10, 430, 131, 23));
        MainWindow->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 816, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        MainWindow->setMenuBar(menuBar);
        QWidget::setTabOrder(header_input, Nonlin_onOff);
        QWidget::setTabOrder(Nonlin_onOff, linear_onOff);
        QWidget::setTabOrder(linear_onOff, initialCondition);
        QWidget::setTabOrder(initialCondition, breaking_beta0);
        QWidget::setTabOrder(breaking_beta0, gravity_input);
        QWidget::setTabOrder(gravity_input, Density);
        QWidget::setTabOrder(Density, timeDuration);
        QWidget::setTabOrder(timeDuration, tStart);
        QWidget::setTabOrder(tStart, workingDir);
        QWidget::setTabOrder(workingDir, OpenDirBrowser);
        QWidget::setTabOrder(OpenDirBrowser, length);
        QWidget::setTabOrder(length, width);
        QWidget::setTabOrder(width, depth);
        QWidget::setTabOrder(depth, nx);
        QWidget::setTabOrder(nx, ny);
        QWidget::setTabOrder(ny, nz);
        QWidget::setTabOrder(nz, sx);
        QWidget::setTabOrder(sx, sy);
        QWidget::setTabOrder(sy, sz);
        QWidget::setTabOrder(sz, alpha);
        QWidget::setTabOrder(alpha, beta);
        QWidget::setTabOrder(beta, gamma);
        QWidget::setTabOrder(gamma, a);
        QWidget::setTabOrder(a, b);
        QWidget::setTabOrder(b, c);
        QWidget::setTabOrder(c, dt);
        QWidget::setTabOrder(dt, waveType);
        QWidget::setTabOrder(waveType, LorP_ComboBox);
        QWidget::setTabOrder(LorP_ComboBox, SF_H);
        QWidget::setTabOrder(SF_H, SF_T);
        QWidget::setTabOrder(SF_T, SF_h);
        QWidget::setTabOrder(SF_h, SF_U);
        QWidget::setTabOrder(SF_U, stokesOrEuler);
        QWidget::setTabOrder(stokesOrEuler, SF_n);
        QWidget::setTabOrder(SF_n, irrFileName);
        QWidget::setTabOrder(irrFileName, xGenStart);
        QWidget::setTabOrder(xGenStart, xGenEnd);
        QWidget::setTabOrder(xGenEnd, yGenStart);
        QWidget::setTabOrder(yGenStart, yGenEnd);
        QWidget::setTabOrder(yGenEnd, rampTime);
        QWidget::setTabOrder(rampTime, xAbsorbStart);
        QWidget::setTabOrder(xAbsorbStart, xAbsorbEnd);
        QWidget::setTabOrder(xAbsorbEnd, yAbsorbStart);
        QWidget::setTabOrder(yAbsorbStart, yAbsorbEnd);
        QWidget::setTabOrder(yAbsorbEnd, pressureDampingOrRelax);
        QWidget::setTabOrder(pressureDampingOrRelax, Hs);
        QWidget::setTabOrder(Hs, Tp);
        QWidget::setTabOrder(Tp, h);
        QWidget::setTabOrder(h, gamma_2);
        QWidget::setTabOrder(gamma_2, maxkh);
        QWidget::setTabOrder(maxkh, seed);
        QWidget::setTabOrder(seed, nOutFiles);
        QWidget::setTabOrder(nOutFiles, storeAscii_onOff);
        QWidget::setTabOrder(storeAscii_onOff, nASCII);
        QWidget::setTabOrder(nASCII, tableWidget);
        QWidget::setTabOrder(tableWidget, pushButton_2);
        QWidget::setTabOrder(pushButton_2, pushButton);
        QWidget::setTabOrder(pushButton, output);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(actionOpen);

        retranslateUi(MainWindow);
        QObject::connect(pushButton, SIGNAL(clicked()), MainWindow, SLOT(close()));

        output->setCurrentIndex(0);
        initialCondition->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "OceanWave3D", 0, QApplication::UnicodeUTF8));
        actionOpen->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("MainWindow", "close", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("MainWindow", "Write input File", 0, QApplication::UnicodeUTF8));
        label_79->setText(QApplication::translate("MainWindow", "Header:", 0, QApplication::UnicodeUTF8));
        label_82->setText(QApplication::translate("MainWindow", "Initial conditions", 0, QApplication::UnicodeUTF8));
        label_87->setText(QApplication::translate("MainWindow", "Breaking filter", 0, QApplication::UnicodeUTF8));
        header_input->setPlainText(QApplication::translate("MainWindow", "Header", 0, QApplication::UnicodeUTF8));
        label_89->setText(QApplication::translate("MainWindow", "Fraction of gravity", 0, QApplication::UnicodeUTF8));
        label_90->setText(QApplication::translate("MainWindow", "Constants", 0, QApplication::UnicodeUTF8));
        label_91->setText(QApplication::translate("MainWindow", "Gravity", 0, QApplication::UnicodeUTF8));
        label_92->setText(QApplication::translate("MainWindow", "Mode", 0, QApplication::UnicodeUTF8));
        linear_onOff->setText(QApplication::translate("MainWindow", "linear", 0, QApplication::UnicodeUTF8));
        Nonlin_onOff->setText(QApplication::translate("MainWindow", "non-linear", 0, QApplication::UnicodeUTF8));
        initialCondition->clear();
        initialCondition->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Still water", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Standing wave", 0, QApplication::UnicodeUTF8)
        );
        label_139->setText(QApplication::translate("MainWindow", "Time", 0, QApplication::UnicodeUTF8));
        label_141->setText(QApplication::translate("MainWindow", "Duration", 0, QApplication::UnicodeUTF8));
        label_145->setText(QApplication::translate("MainWindow", "Start time", 0, QApplication::UnicodeUTF8));
        label_158->setText(QApplication::translate("MainWindow", "Density", 0, QApplication::UnicodeUTF8));
        label_159->setText(QApplication::translate("MainWindow", "m/s^2", 0, QApplication::UnicodeUTF8));
        label_160->setText(QApplication::translate("MainWindow", "kg/m^3", 0, QApplication::UnicodeUTF8));
        label_161->setText(QApplication::translate("MainWindow", "s", 0, QApplication::UnicodeUTF8));
        label_162->setText(QApplication::translate("MainWindow", "s", 0, QApplication::UnicodeUTF8));
        label_165->setText(QApplication::translate("MainWindow", "Working directory ", 0, QApplication::UnicodeUTF8));
        OpenDirBrowser->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        workingDir->setText(QString());
        output->setTabText(output->indexOf(tab_11), QApplication::translate("MainWindow", "General", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("MainWindow", "nx", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("MainWindow", "Stretching", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "depth [m]", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_WHATSTHIS
        label_2->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        label_2->setText(QApplication::translate("MainWindow", "length [m]", 0, QApplication::UnicodeUTF8));
        sy->setText(QApplication::translate("MainWindow", "y", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("MainWindow", "nz", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("MainWindow", "Resolution", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("MainWindow", "ny", 0, QApplication::UnicodeUTF8));
        sz->setText(QApplication::translate("MainWindow", "z", 0, QApplication::UnicodeUTF8));
        sx->setText(QApplication::translate("MainWindow", "x", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("MainWindow", "Domain", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "width [m]", 0, QApplication::UnicodeUTF8));
        label_88->setText(QApplication::translate("MainWindow", "Geometry", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(tab), QApplication::translate("MainWindow", "Geometry", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("MainWindow", "Discretization", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("MainWindow", "Preconditioner", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("MainWindow", "Finite Difference", 0, QApplication::UnicodeUTF8));
        label_27->setText(QApplication::translate("MainWindow", "b", 0, QApplication::UnicodeUTF8));
        label_28->setText(QApplication::translate("MainWindow", "gamma", 0, QApplication::UnicodeUTF8));
        label_29->setText(QApplication::translate("MainWindow", "beta", 0, QApplication::UnicodeUTF8));
        label_30->setText(QApplication::translate("MainWindow", "a", 0, QApplication::UnicodeUTF8));
        label_31->setText(QApplication::translate("MainWindow", "c", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_WHATSTHIS
        label_32->setWhatsThis(QString());
#endif // QT_NO_WHATSTHIS
        label_32->setText(QApplication::translate("MainWindow", "alpha", 0, QApplication::UnicodeUTF8));
        label_137->setText(QApplication::translate("MainWindow", "Time", 0, QApplication::UnicodeUTF8));
        label_138->setText(QApplication::translate("MainWindow", "dt", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(tab_2), QApplication::translate("MainWindow", "Discretization", 0, QApplication::UnicodeUTF8));
        label_52->setText(QApplication::translate("MainWindow", "Wave generation and absorption", 0, QApplication::UnicodeUTF8));
        label_55->setText(QApplication::translate("MainWindow", "Wave generation zone", 0, QApplication::UnicodeUTF8));
        label_126->setText(QApplication::translate("MainWindow", "x [m]", 0, QApplication::UnicodeUTF8));
        label_127->setText(QApplication::translate("MainWindow", "y [m]", 0, QApplication::UnicodeUTF8));
        label_129->setText(QApplication::translate("MainWindow", "Begin", 0, QApplication::UnicodeUTF8));
        label_130->setText(QApplication::translate("MainWindow", "End", 0, QApplication::UnicodeUTF8));
        label_132->setText(QApplication::translate("MainWindow", "x [m]", 0, QApplication::UnicodeUTF8));
        label_133->setText(QApplication::translate("MainWindow", "End", 0, QApplication::UnicodeUTF8));
        label_134->setText(QApplication::translate("MainWindow", "Wave absoption zone", 0, QApplication::UnicodeUTF8));
        label_135->setText(QApplication::translate("MainWindow", "y[m]", 0, QApplication::UnicodeUTF8));
        label_136->setText(QApplication::translate("MainWindow", "Begin", 0, QApplication::UnicodeUTF8));
        SF_TLabel->setText(QApplication::translate("MainWindow", "Wave period:", 0, QApplication::UnicodeUTF8));
        label_142->setText(QApplication::translate("MainWindow", "Water depth", 0, QApplication::UnicodeUTF8));
        label_140->setText(QApplication::translate("MainWindow", "Wave height", 0, QApplication::UnicodeUTF8));
        label_143->setText(QApplication::translate("MainWindow", "Current velocity", 0, QApplication::UnicodeUTF8));
        label_144->setText(QApplication::translate("MainWindow", "m", 0, QApplication::UnicodeUTF8));
        SF_TL_unit->setText(QApplication::translate("MainWindow", "s", 0, QApplication::UnicodeUTF8));
        label_146->setText(QApplication::translate("MainWindow", "m", 0, QApplication::UnicodeUTF8));
        label_147->setText(QApplication::translate("MainWindow", "m/s", 0, QApplication::UnicodeUTF8));
        label_148->setText(QApplication::translate("MainWindow", "Current type", 0, QApplication::UnicodeUTF8));
        stokesOrEuler->clear();
        stokesOrEuler->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Euler", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Stokes", 0, QApplication::UnicodeUTF8)
        );
        label_149->setText(QApplication::translate("MainWindow", "Order", 0, QApplication::UnicodeUTF8));
        irrFileName->setText(QApplication::translate("MainWindow", "DummyName", 0, QApplication::UnicodeUTF8));
        irrFileName_label->setText(QApplication::translate("MainWindow", "File name", 0, QApplication::UnicodeUTF8));
        LorP_ComboBox->clear();
        LorP_ComboBox->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Wave length", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "wave period", 0, QApplication::UnicodeUTF8)
        );
        pressureDampingOrRelax->clear();
        pressureDampingOrRelax->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Pressure damping", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Relaxation zone", 0, QApplication::UnicodeUTF8)
        );
        label_150->setText(QApplication::translate("MainWindow", "Tp", 0, QApplication::UnicodeUTF8));
        label_151->setText(QApplication::translate("MainWindow", "Water depth", 0, QApplication::UnicodeUTF8));
        label_152->setText(QApplication::translate("MainWindow", "Hs", 0, QApplication::UnicodeUTF8));
        label_153->setText(QApplication::translate("MainWindow", "gamma", 0, QApplication::UnicodeUTF8));
        label_154->setText(QApplication::translate("MainWindow", "m", 0, QApplication::UnicodeUTF8));
        label_155->setText(QApplication::translate("MainWindow", "s", 0, QApplication::UnicodeUTF8));
        label_156->setText(QApplication::translate("MainWindow", "m", 0, QApplication::UnicodeUTF8));
        label_157->setText(QApplication::translate("MainWindow", "m/s", 0, QApplication::UnicodeUTF8));
        label_166->setText(QApplication::translate("MainWindow", "max(kh)", 0, QApplication::UnicodeUTF8));
        label_167->setText(QApplication::translate("MainWindow", "seed", 0, QApplication::UnicodeUTF8));
        label_163->setText(QApplication::translate("MainWindow", "Ramp time [s]", 0, QApplication::UnicodeUTF8));
        label_164->setText(QApplication::translate("MainWindow", "type", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(tab_5), QApplication::translate("MainWindow", "Wave Generation", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Number of output files", 0, QApplication::UnicodeUTF8));
        storeAscii_onOff->setText(QApplication::translate("MainWindow", "Store ASCII files", 0, QApplication::UnicodeUTF8));
        ACCII_label->setText(QApplication::translate("MainWindow", "ASCII files are stored for every:", 0, QApplication::UnicodeUTF8));
        ASCII_label2->setText(QApplication::translate("MainWindow", "time steps.", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(tab_9), QApplication::translate("MainWindow", "Ouput", 0, QApplication::UnicodeUTF8));
        selectPPfile->setText(QApplication::translate("MainWindow", "Select File", 0, QApplication::UnicodeUTF8));
        convert->setText(QApplication::translate("MainWindow", "Convert", 0, QApplication::UnicodeUTF8));
        selectGnuplotFile->setText(QApplication::translate("MainWindow", "Select File", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("MainWindow", "GNUPLOT", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("MainWindow", "Convert", 0, QApplication::UnicodeUTF8));
        plot->setText(QApplication::translate("MainWindow", "plot", 0, QApplication::UnicodeUTF8));
        read_bottom->setText(QApplication::translate("MainWindow", "read", 0, QApplication::UnicodeUTF8));
        SelectOutput->clear();
        SelectOutput->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Select output", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Force", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "matlab", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Raw", 0, QApplication::UnicodeUTF8)
        );
        label_16->setText(QApplication::translate("MainWindow", "x0 [m]", 0, QApplication::UnicodeUTF8));
        label_17->setText(QApplication::translate("MainWindow", "D [m]", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("MainWindow", "Cd", 0, QApplication::UnicodeUTF8));
        label_19->setText(QApplication::translate("MainWindow", "Cm", 0, QApplication::UnicodeUTF8));
        convertStatus->setText(QApplication::translate("MainWindow", "Converting...", 0, QApplication::UnicodeUTF8));
        convertWarning->setPlainText(QApplication::translate("MainWindow", "WARNING: Do not convert files before computation has finished", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(tab_3), QApplication::translate("MainWindow", "Post-processing", 0, QApplication::UnicodeUTF8));
        about_combobox->clear();
        about_combobox->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "OceanWave3D", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "GUI", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Publications", 0, QApplication::UnicodeUTF8)
        );
        aboutText_OCW3D->setHtml(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Ubuntu'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">OceanWave3D - a costal engineering tool for simulation of nonlinear free</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">surface waves. </p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Copyright (C) 2009 Allan P. Engsig-Karup, DTU Compute, </p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Technical University of Denmark</p>\n"
""
                        "<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">This OceanWave3D program comes with ABSOLUTELY NO WARRANTY. This is free</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">software and you are welcome to redistribute it under the conditions</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">of the GNU General Public License version 3.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Permissi"
                        "on is hereby granted, free of charge, to any person obtaining a copy</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">of this software and associated documentation files (the &quot;Software&quot;), to deal</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">in the Software without restriction, including without limitation the rights</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">to use, copy, modify, merge, publish, distribute, sublicense, and/or sell</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">copies of the Software, and to permit persons to whom the Software is</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">fu"
                        "rnished to do so, subject to the following conditions:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">The above copyright notice and this permission notice shall be included in</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">all copies or substantial portions of the Software.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">THE SOFTWARE IS PROVIDED &quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR</p>\n"
"<p style=\" margin-top:0px; margin-bottom"
                        ":0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; ma"
                        "rgin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">THE SOFTWARE.</p></body></html>", 0, QApplication::UnicodeUTF8));
        aboutText_OCW3dGUI->setHtml(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Ubuntu'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">The Graphical User Interface (GUI) for OceanWave3D (C) is developed and maintaied by Bo Terp Paulsen, bo.paulsen@deltares.nl. The software has inherrited the GNU General Public Lisence version 3 of OceanWave3D, and is hence free to use and redistribute. </p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">THE SOFTWARE IS PROVIDED &"
                        "quot;AS IS&quot;, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">OUT OF OR IN CONNECTION"
                        " WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">THE SOFTWARE.</p></body></html>", 0, QApplication::UnicodeUTF8));
        aboutText_OCW3D_publications->setHtml(QApplication::translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:'Ubuntu'; font-size:9pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Please acknowledge this software referred to as &quot;OceanWave3D&quot; in your work and any of the publications below by using the following references:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  @ARTICLE{EngsigKarupEtAl08,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-blo"
                        "ck-indent:0; text-indent:0px;\">  AUTHOR    = &quot;Engsig-Karup, A.P. and Bingham, H.B. and Lindberg, O.&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  TITLE     = &quot;An efficient flexible-order model for {3D} nonlinear water waves&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  YEAR      = &quot;2009&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  JOURNAL   = &quot;Journal of Computational Physics&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  VOLUME    = &quot;228&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  PAGES     = &quot;2100-2118&quot;</p>\n"
"<p style"
                        "=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">and this chapter with general details on the properties of the model and its parallel implementation (note that this Fortran version is not yet parallel)</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  @inbook{EngsigKarupEtAl2013,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  title     = "
                        "&quot;Fast hydrodynamics on heterogenous many-core hardware&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  publisher = &quot;Taylor &amp; Francis&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  author    = &quot;Engsig-Karup, {Allan Peter} and Glimberg, {Stefan Lemvig} and Nielsen, {Allan S.} and Ole Lindberg&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  note      = &quot;2013; 11&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  year      = &quot;2013&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  editor    = &quot;Rapha\\'el Couturier&quot;,</p>\n"
"<p style=\" mar"
                        "gin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  isbn      = &quot;978-1-4665-7162-4&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  pages     = &quot;251--294&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  booktitle = &quot;Designing Scientific Applications on GPUs&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">  }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">OceanWave3D is distributed under the GNU General Public License (See file LICENSE) and a"
                        " base code was developed at DTU Mechanics 2006-2008 by Allan P. Engsig-Karup and entirely rewritten at DTU Informatics 2008-2011 (current version) by Allan P. Engsig-Karup with contributions as mentioned below section 4 of this file.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">A GPU-accelerated massively parallel version of OceanWave3D that can execute on heterogenous archictures using the CUDA programming model has been developed in collaboration between Morten Gorm Madsen (first proof-of-concept version), Stefan Lemvig Glimberg (advanced version within GPULAB library) and Allan P. Engsig-Karup. In this version of the code the multigrid preconditioned defect correction method (PDC) appeared for the first time and it accessbile in this fortran version of the code as we"
                        "ll. This work can be cited by the following reference</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   @ARTICLE{EngsigKarupEtAl2011,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   AUTHOR = &quot;Engsig-Karup, A. P. and Madsen, M. G. and Glimberg, S. L.&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   TITLE = &quot;A massively parallel {GPU}-accelerated model for analysis of fully nonlinear free surface waves&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   JOURNAL = &quot;International Journal of Numerical Meth"
                        "ods in Fluids&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   YEAR = &quot;2011&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   DOI = &quot;10.1002/fld.2675&quot;</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">The properties of the PDC method is analyzed in great detail in</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:"
                        "0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   @ARTICLE{EngsigKarupEtAl2013,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   AUTHOR = &quot;Engsig-Karup, A.P.&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   TITLE = &quot;Analysis of Efficient Preconditioned Defect Correction Methods for Nonlinear Water Waves&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   JOURNAL = &quot;International Journal of Numerical Methods in Fluids&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   YEAR = &quot;2013&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text"
                        "-indent:0px;\">   }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">The splitting technique that can be used for efficient wave-structure interactions are described here</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   @article{DucrozetEtAl2013,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   title     = &quot;A non-linear wave decomposition model for efficient wave\342\200\232\303\204\303\254structure interaction. Part A: Formulation, validations and"
                        " analysis&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   publisher = &quot;Academic Press&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   author    = &quot;Guillaume Ducrozet and Engsig-Karup, {Allan Peter} and Bingham, {Harry B.} and Pierre Ferrant&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   year      = &quot;2014&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   volume    = &quot;257&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   pages     = &quot;863--883&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-bl"
                        "ock-indent:0; text-indent:0px;\">   journal   = &quot;Journal of Computational Physics&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   issn      = &quot;0021-9991&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">A comparison between OceanWave3D and HOS is described here</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-inden"
                        "t:0px;\">   @article{DucrozetEtAl2012,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   title     = &quot;A comparative study of two fast nonlinear free-surface water wave models&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   publisher = &quot;John/Wiley &amp; Sons Ltd.&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   author    = &quot;Guillaume Ducrozet and Bingham, {Harry B.} and Engsig-Karup, {Allan Peter} and Felicien Bonnefoy and Pierre Ferrant&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   year      = &quot;2012&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   volum"
                        "e    = &quot;69&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   pages     = &quot;1818--1834&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   journal   = &quot;International Journal for Numerical Methods in Fluids&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   issn      = &quot;0271-2091&quot;,</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">   }</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>", 0, QApplication::UnicodeUTF8));
        output->setTabText(output->indexOf(About), QApplication::translate("MainWindow", "About", 0, QApplication::UnicodeUTF8));
        run->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
