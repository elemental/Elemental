/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_DISPLAYWINDOW_DECL_HPP
#define GRAPHICS_DISPLAYWINDOW_DECL_HPP

#include <QPushButton>
#include <QScrollArea>
#include <QWidget>

#include "elemental-lite.hpp"
#include "elemental/graphics/display_widget_decl.hpp"

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
    void UseLocalScale();
    void UseGlobalScale();
};

} // namespace elem

#endif // ifndef GRAPHICS_DISPLAYWINDOW_DECL_HPP
