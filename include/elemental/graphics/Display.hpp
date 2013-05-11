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

template<typename T>
DisplayWidget<T>::DisplayWidget( QWidget* parent )
: QWidget(parent)
{ }

template<typename T>
DisplayWidget<T>::~DisplayWidget()
{ }

template<typename T>
void DisplayWidget<T>::paintEvent( QPaintEvent* event )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::paintEvent");
#endif
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap );
}

template<typename T>
void DisplayWidget<T>::DisplayReal( const Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayReal");
#endif
    typedef BASE(T) R;
    const int m = A.Height();
    const int n = A.Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, m );
    const int nPix = std::max( 500, n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Compute the range of the real values in A
    R minVal=0, maxVal=0;
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A.GetRealPart( 0, 0 );
        for( int j=0; j<n; ++j )
        {
            for( int i=0; i<m; ++i )
            {
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
            QRgb color = ColorMap( A.GetRealPart(i,j), minVal, maxVal );
            painter.setPen( color );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();
}

template<typename T>
void DisplayWidget<T>::DisplayImag( const Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayImag");
#endif
    typedef BASE(T) R;
    const int m = A.Height();
    const int n = A.Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, m );
    const int nPix = std::max( 500, n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Compute the range of the real values in A
    R minVal=0, maxVal=0;
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A.GetImagPart( 0, 0 );
        for( int j=0; j<n; ++j )
        {
            for( int i=0; i<m; ++i )
            {
                minVal = std::min( minVal, A.GetImagPart(i,j) );
                maxVal = std::max( maxVal, A.GetImagPart(i,j) );
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
            QRgb color = ColorMap( A.GetImagPart(i,j), minVal, maxVal );
            painter.setPen( color );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();
}

template<typename T>
DisplayWindow<T>::DisplayWindow( QString title, QWidget* parent )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::DisplayWindow");
#endif
    QHBoxLayout* layout = new QHBoxLayout( this );

    // Display the real data
    realDisplay = new DisplayWidget<T>();
    realScroll = new QScrollArea();
    realScroll->setWidget( realDisplay );
    layout->addWidget( realScroll );

    // If the data is complex, also show the imaginary part
    if( IsComplex<T>::val )
    {
        imagDisplay = new DisplayWidget<T>();
        imagScroll = new QScrollArea();
        imagScroll->setWidget( imagDisplay );
        layout->addWidget( imagScroll );
    }

    setLayout( layout );
    setWindowTitle( title );
    show();

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

template<typename T>
DisplayWindow<T>::DisplayWindow
( const Matrix<T>& A, QString title, QWidget* parent )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::DisplayWindow");
#endif
    QHBoxLayout* layout = new QHBoxLayout( this );

    // Display the real data
    realDisplay = new DisplayWidget<T>();
    realScroll = new QScrollArea();
    realScroll->setWidget( realDisplay );
    layout->addWidget( realScroll );

    // If the data is complex, also show the imaginary part
    if( IsComplex<T>::val )
    {
        imagDisplay = new DisplayWidget<T>();
        imagScroll = new QScrollArea();
        imagScroll->setWidget( imagDisplay );
        layout->addWidget( imagScroll );
    }

    setLayout( layout );
    Display( A, title );
    show();

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

template<typename T>
DisplayWindow<T>::~DisplayWindow()
{ }

template<typename T>
void DisplayWindow<T>::Display( const Matrix<T>& A, QString title )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::Display");
#endif
    setWindowTitle( title );
    realDisplay->DisplayReal( A );
    if( IsComplex<T>::val )
        imagDisplay->DisplayImag( A );
}

template<typename T>
void Display( const Matrix<T>& A, std::string title="" )
{
#ifndef RELEASE
    CallStackEntry entry("Display");
#endif
    QString qTitle = QString::fromStdString( title );
    DisplayWindow<T>* displayWindow = new DisplayWindow<T>( A, qTitle );
    RegisterWindow( displayWindow );
}

template<typename T,Distribution U,Distribution V>
void Display( const DistMatrix<T,U,V>& A, std::string title="" )
{
#ifndef RELEASE
    CallStackEntry entry("Display");
#endif
    // TODO: Avoid giving every process a full copy
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    if( A.Grid().Rank() == 0 )
        Display( A_STAR_STAR.Matrix(), title );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_DISPLAY_HPP
