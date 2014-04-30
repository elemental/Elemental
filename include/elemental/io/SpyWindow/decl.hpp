/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SPYWINDOW_DECL_HPP
#define ELEM_SPYWINDOW_DECL_HPP
#ifdef ELEM_HAVE_QT5

#include <QScrollArea>
#include <QWidget>

#include "elemental-lite.hpp"
#include "elemental/io/SpyWidget/decl.hpp"

namespace elem {

class SpyWindow : public QWidget
{
public:
    SpyWindow( QWidget* parent=0 );    
    ~SpyWindow();

    void Spy
    ( const Matrix<Int>* A, 
      QString title=QString("Default title") );

private:
    QScrollArea *scroll_;
    SpyWidget *spy_;
    const Matrix<Int> *matrix_;
};

} // namespace elem

#endif // ifdef ELEM_HAVE_QT5
#endif // ifndef ELEM_SPYWINDOW_DECL_HPP
