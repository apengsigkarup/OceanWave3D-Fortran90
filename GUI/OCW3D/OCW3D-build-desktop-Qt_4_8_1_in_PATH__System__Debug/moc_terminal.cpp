/****************************************************************************
** Meta object code from reading C++ file 'terminal.h'
**
** Created: Tue Apr 14 13:33:34 2015
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../OCW3D/terminal.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'terminal.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_Terminal[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x05,

 // slots: signature, parameters, type, tag, flags
      24,    9,   19,    9, 0x0a,
      32,    9,   19,    9, 0x0a,
      47,    9,    9,    9, 0x09,

       0        // eod
};

static const char qt_meta_stringdata_Terminal[] = {
    "Terminal\0\0exited()\0bool\0start()\0"
    "tryTerminate()\0termProcessExited()\0"
};

void Terminal::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        Terminal *_t = static_cast<Terminal *>(_o);
        switch (_id) {
        case 0: _t->exited(); break;
        case 1: { bool _r = _t->start();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 2: { bool _r = _t->tryTerminate();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 3: _t->termProcessExited(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData Terminal::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject Terminal::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_Terminal,
      qt_meta_data_Terminal, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &Terminal::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *Terminal::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *Terminal::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_Terminal))
        return static_cast<void*>(const_cast< Terminal*>(this));
    return QWidget::qt_metacast(_clname);
}

int Terminal::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void Terminal::exited()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
