/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_SPYWINDOW_DECL_HPP
#define ELEM_IO_SPYWINDOW_DECL_HPP

#include <QScrollArea>
#include <QWidget>

#include "elemental-lite.hpp"
#include "elemental/io/spy_widget_decl.hpp"

namespace elem {

class SpyWindow : public QWidget
{
public:
    SpyWindow( QWidget* parent=0 );    
    ~SpyWindow();

    void Spy
    ( const Matrix<int>* A, 
      QString title=QString("Default title") );

private:
    QScrollArea *scroll_;
    SpyWidget *spy_;
    const Matrix<int> *matrix_;
};

} // namespace elem

#endif // ifndef ELEM_IO_SPYWINDOW_DECL_HPP
