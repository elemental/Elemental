/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPYWINDOW_HPP
#define EL_SPYWINDOW_HPP
#ifdef EL_HAVE_QT5

#include <QScrollArea>
#include <QWidget>

#include "El.hpp"

#include "El/io/SpyWidget.hpp"

namespace El {

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

} // namespace El

#endif // ifdef EL_HAVE_QT5
#endif // ifndef EL_SPYWINDOW_HPP
