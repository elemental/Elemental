/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IO_DISPLAYWINDOW_DECL_HPP
#define EL_IO_DISPLAYWINDOW_DECL_HPP

// Not currently guarded due to CMake MOC handling requiring extra flags
//#ifdef EL_HAVE_QT5

#include <QPushButton>
#include <QScrollArea>
#include <QWidget>

#include "El/config.h"
#ifdef EL_HAVE_F90_INTERFACE
# include "El/FCMangle.h"
#endif
#include "El/core.hpp"
#include "El/io.hpp"

namespace El {

// Unfortunately Q_OBJECT does not support templates...
class DisplayWindow : public QWidget
{
    Q_OBJECT

public:
    DisplayWindow( QWidget* parent=0 );    
    ~DisplayWindow();

    void Display
    ( const Matrix<double>* A, 
      QString title=QString("Default title") );
    void Display
    ( const Matrix<double>* A, 
      double minVal, double maxVal,
      QString title=QString("Default title") );

private:
    QScrollArea *scroll_;
    DisplayWidget<double> *display_;
    const Matrix<double> *matrix_;

public slots:
    void Save();
    void SetScale( bool global );
};

} // namespace El

//#endif // ifdef EL_HAVE_QT5

#endif // ifndef EL_IO_DISPLAYWINDOW_DECL_HPP
