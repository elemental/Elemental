/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_DISPLAY_HPP
#define GRAPHICS_DISPLAY_HPP

#ifdef HAVE_QT5

#include "elemental/graphics/ColorMap.hpp"

namespace elem {

template<typename R>
DisplayWidget<R>::DisplayWidget( QWidget* parent )
: QWidget(parent)
{ }

template<typename R>
DisplayWidget<R>::~DisplayWidget()
{ }

template<typename R>
void DisplayWidget<R>::paintEvent( QPaintEvent* event )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::paintEvent");
#endif
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap );
}

template<typename R>
void DisplayWidget<R>::Display( const Matrix<R>& A )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::Display");
#endif
    const int m = A.Height();
    const int n = A.Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, m );
    const int nPix = std::max( 500, n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Compute the range of the values in A
    BASE(R) minVal=0, maxVal=0;
    if( m != 0 && n != 0 )
    {
        //minVal = maxVal = A.Get( 0, 0 );
        minVal = maxVal = A.GetRealPart( 0, 0 );
        for( int j=0; j<n; ++j )
        {
            for( int i=0; i<m; ++i )
            {
                //minVal = std::min( minVal, A.Get(i,j) );
                //maxVal = std::max( maxVal, A.Get(i,j) );
                minVal = std::min( minVal, A.GetRealPart(i,j) );
                maxVal = std::max( maxVal, A.GetRealPart(i,j) );
            }
        }
    }

    // Paint the matrix
    QPainter painter( &pixmap );
    painter.initFrom( this );
    for( int jPix=0; jPix<nPix; ++jPix ) 
    {
        const int j = nRatio*jPix;
        for( int iPix=0; iPix<mPix; ++iPix )
        {
            const int i = mRatio*iPix;
            //QRgb color = ColorMap( A.Get(i,j), minVal, maxVal );
            QRgb color = ColorMap( A.GetRealPart(i,j), minVal, maxVal );
            painter.setPen( color );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();
}

template<typename T>
DisplayWindow<T>::DisplayWindow( QWidget* parent )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::DisplayWindow");
#endif
    display = new DisplayWidget<T>();
    scroll = new QScrollArea();
    scroll->setWidget( display );

    QVBoxLayout* layout = new QVBoxLayout( this );
    layout->addWidget( scroll );
    setLayout( layout );

    OpenedQtWindow();
}

template<typename T>
DisplayWindow<T>::DisplayWindow( const Matrix<T>& A, QWidget* parent )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::DisplayWindow");
#endif
    display = new DisplayWidget<T>();
    scroll = new QScrollArea();
    scroll->setWidget( display );

    QVBoxLayout* layout = new QVBoxLayout( this );
    layout->addWidget( scroll );
    setLayout( layout );

    show();
    display->Display( A );

    OpenedQtWindow();
}

template<typename T>
DisplayWindow<T>::~DisplayWindow()
{ }

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_DISPLAY_HPP
