/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Thu Apr 23 14:57:53 2015
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../OCW3D/mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      49,   38,   11,   11, 0x08,
      75,   11,   11,   11, 0x08,
     105,   98,   11,   11, 0x08,
     141,  133,   11,   11, 0x08,
     158,   11,   11,   11, 0x08,
     175,   11,   11,   11, 0x08,
     195,   11,   11,   11, 0x08,
     210,   11,   11,   11, 0x08,
     216,   11,   11,   11, 0x08,
     226,   11,   11,   11, 0x08,
     246,   11,   11,   11, 0x08,
     265,   11,   11,   11, 0x08,
     279,   11,   11,   11, 0x08,
     291,   11,   11,   11, 0x08,
     312,   11,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0on_pushButton_2_clicked()\0"
    "waveTheory\0on_waveTheoryChanged(int)\0"
    "on_waveTheoryChanged()\0nFiles\0"
    "on_outputWidgetChanged(int)\0checked\0"
    "storeASCII(bool)\0openFileDialog()\0"
    "openWorkDirDialog()\0selectPPfile()\0"
    "run()\0gnuplot()\0readKinematicFile()\0"
    "about_changed(int)\0saveAsAscii()\0"
    "saveAsMat()\0convertTo_setup(int)\0"
    "convertTo()\0"
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        MainWindow *_t = static_cast<MainWindow *>(_o);
        switch (_id) {
        case 0: _t->on_pushButton_2_clicked(); break;
        case 1: _t->on_waveTheoryChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->on_waveTheoryChanged(); break;
        case 3: _t->on_outputWidgetChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: _t->storeASCII((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 5: _t->openFileDialog(); break;
        case 6: _t->openWorkDirDialog(); break;
        case 7: _t->selectPPfile(); break;
        case 8: _t->run(); break;
        case 9: _t->gnuplot(); break;
        case 10: _t->readKinematicFile(); break;
        case 11: _t->about_changed((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: _t->saveAsAscii(); break;
        case 13: _t->saveAsMat(); break;
        case 14: _t->convertTo_setup((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: _t->convertTo(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData MainWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
