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

#include <QBoxLayout>
#include <QPainter>
#include <QPixmap>
#include <QScrollArea>
#include <QStylePainter>
#include <QWidget>

namespace elem {

template<typename T>
QRgb ColorMap( T value, T minVal, T maxVal );

template<typename T>
class DisplayWidget : public QWidget
{
public:
    DisplayWidget( QWidget* parent=0 );
    ~DisplayWidget();
    // TODO: Generalize to function which displays f(A), where f is functor
    void DisplayReal( const Matrix<T>& A );
    void DisplayImag( const Matrix<T>& A );
    // TODO: Add colorbar
protected:
    void paintEvent( QPaintEvent* event );
private:
    QPixmap pixmap;
};

// This exists so that we can push all windows 
// (with potentially different matrix datatype) to the same stack
class Window : public QWidget
{ };

template<typename T> 
class DisplayWindow : public Window
{
public:
    explicit DisplayWindow
    ( QString title=QString("Default title"), QWidget* parent=0 );    
    explicit DisplayWindow
    ( const Matrix<T>& A, 
      QString title=QString("Default title"), QWidget* parent=0 );    
    ~DisplayWindow();

    void Display
    ( const Matrix<T>& A, QString title=QString("Default title") );
private:
    QScrollArea *realScroll, *imagScroll;
    DisplayWidget<T> *realDisplay, *imagDisplay;
};

// For letting Elemental handle deleting windows in 'Finalize'
void RegisterWindow( Window* window );

// When Elemental is finalized, if no window was opened, then it must call 
// app.exit() instead
void OpenedWindow();

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_DECL_HPP
