/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_DISPLAYWINDOW_DECL_HPP
#define ELEM_IO_DISPLAYWINDOW_DECL_HPP

// Not currently guarded due to CMake MOC handling requiring extra flags
//#ifdef ELEM_HAVE_QT5

#include <QPushButton>
#include <QScrollArea>
#include <QWidget>

#include "elemental/include-paths.hpp"
#include "elemental/config.h"
#ifdef ELEM_HAVE_F90_INTERFACE
# include "elemental/FCMangle.h"
#endif
#include "elemental/core.hpp"
#include "elemental/blas-like/decl.hpp"
#include "elemental/lapack-like/decl.hpp"
#include "elemental/convex/decl.hpp"
#include "elemental/io/DisplayWidget/decl.hpp"

namespace elem {

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

} // namespace elem

//#endif // ifdef ELEM_HAVE_QT5

#endif // ifndef ELEM_IO_DISPLAYWINDOW_DECL_HPP
