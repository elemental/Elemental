/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_DECL_HPP
#define GRAPHICS_DECL_HPP

#ifdef HAVE_QT5

#include <QPainter>
#include <QPixmap>
#include <QScrollArea>
#include <QStylePainter>
#include <QVBoxLayout>
#include <QWidget>

namespace elem {

// When Elemental is finalized, if no window was opened, then it must call 
// app.exit() instead
void OpenedQtWindow();

template<typename R>
QRgb ColorMap( R value, R minVal, R maxVal );

template<typename R>
class DisplayWidget : public QWidget
{
public:
    DisplayWidget( QWidget* parent=0 );
    ~DisplayWidget();
    void Display( const Matrix<R>& A );
protected:
    void paintEvent( QPaintEvent* event );
private:
    QPixmap pixmap;
};

/*
template<typename R>
class DisplayWidget<Complex<R> > : public QWidget
{
public:
    DisplayWidget( QWidget* parent=0 );
    ~DisplayWidget();
    void Display( const Matrix<Complex<R> >& A );
protected:
    void paintEvent( QPaintEvent* event );
private:
    QPixmap pixmap;
};
*/

template<typename F> 
class DisplayWindow : public QWidget
{
public:
    explicit DisplayWindow( QWidget* parent=0 );    
    explicit DisplayWindow( const Matrix<F>& A, QWidget* parent=0 );    
    ~DisplayWindow();
    void Display( const Matrix<F>& A ) const;
private:
    QScrollArea* scroll;
    DisplayWidget<F>* display;    
};

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_DECL_HPP
