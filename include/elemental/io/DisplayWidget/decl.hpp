/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DISPLAYWIDGET_DECL_HPP
#define ELEM_DISPLAYWIDGET_DECL_HPP
#ifdef ELEM_HAVE_QT5

#include <QWidget>

namespace elem {

template<typename T>
class DisplayWidget : public QWidget
{
public:
    DisplayWidget( QWidget* parent=0 );
    ~DisplayWidget();
    // TODO: Generalize to function which displays f(A), where f is functor
    void DisplayReal( const Matrix<T>* A );
    void DisplayReal( const Matrix<T>* A, Base<T> minVal, Base<T> maxVal );
    void DisplayImag( const Matrix<T>* A );
    void DisplayImag( const Matrix<T>* A, Base<T> minVal, Base<T> maxVal );
    // TODO: Add colorbar

    void SavePng( std::string basename ) const;

protected:
    void paintEvent( QPaintEvent* event );

private:
    QPixmap pixmap_;
};

} // namespace elem

#endif // ifdef ELEM_HAVE_QT5
#endif // ifndef ELEM_DISPLAYWIDGET_DECL_HPP
