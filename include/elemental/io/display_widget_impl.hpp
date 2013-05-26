/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef IO_DISPLAYWIDGET_IMPL_HPP
#define IO_DISPLAYWIDGET_IMPL_HPP

#ifdef HAVE_QT5

#include <QFile>
#include <QPainter>
#include <QPixmap>
#include <QStylePainter>

#include "elemental/io/ColorMap.hpp"

#include "elemental/io/display_widget_decl.hpp"

namespace elem {

template<typename T>
inline
DisplayWidget<T>::DisplayWidget( QWidget* parent )
: QWidget(parent)
{ }

template<typename T>
inline
DisplayWidget<T>::~DisplayWidget()
{ }

template<typename T>
inline void 
DisplayWidget<T>::paintEvent( QPaintEvent* event )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::paintEvent");
#endif
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap_ );
}

template<typename T>
inline void 
DisplayWidget<T>::DisplayReal( const Matrix<T>* A )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayReal");
#endif
    typedef BASE(T) R;
    const int m = A->Height();
    const int n = A->Width();

    // Compute the range of the real values in A
    R minVal=0, maxVal=0;
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A->GetRealPart( 0, 0 );
        for( int j=0; j<n; ++j )
        {
            for( int i=0; i<m; ++i )
            {
                minVal = std::min( minVal, A->GetRealPart(i,j) );
                maxVal = std::max( maxVal, A->GetRealPart(i,j) );
            }
        }
    }

    DisplayReal( A, minVal, maxVal );
}

template<typename T>
inline void 
DisplayWidget<T>::DisplayReal
( const Matrix<T>* A, BASE(T) minVal, BASE(T) maxVal )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayReal");
#endif
    typedef BASE(T) R;
    const int m = A->Height();
    const int n = A->Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, 2*m );
    const int nPix = std::max( 500, 2*n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( int jPix=0; jPix<nPix; ++jPix ) 
    {
        const int j = nRatio*jPix;
        for( int iPix=0; iPix<mPix; ++iPix )
        {
            const int i = mRatio*iPix;
            QRgb color = ColorMap( A->GetRealPart(i,j), minVal, maxVal );
            painter.setPen( color );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();

    // Keep track of the extrema to allow for consistent visualization
    UpdateMinRealWindowVal( minVal );
    UpdateMaxRealWindowVal( maxVal );
}

template<typename T>
inline void 
DisplayWidget<T>::DisplayImag( const Matrix<T>* A )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayImag");
#endif
    typedef BASE(T) R;
    const int m = A->Height();
    const int n = A->Width();

    // Compute the range of the real values in A
    R minVal=0, maxVal=0;
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A->GetImagPart( 0, 0 );
        for( int j=0; j<n; ++j )
        {
            for( int i=0; i<m; ++i )
            {
                minVal = std::min( minVal, A->GetImagPart(i,j) );
                maxVal = std::max( maxVal, A->GetImagPart(i,j) );
            }
        }
    }

    DisplayImag( A, minVal, maxVal );
}

template<typename T>
inline void 
DisplayWidget<T>::DisplayImag
( const Matrix<T>* A, BASE(T) minVal, BASE(T) maxVal )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::DisplayImag");
#endif
    typedef BASE(T) R;
    const int m = A->Height();
    const int n = A->Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, 2*m );
    const int nPix = std::max( 500, 2*n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( int jPix=0; jPix<nPix; ++jPix ) 
    {
        const int j = nRatio*jPix;
        for( int iPix=0; iPix<mPix; ++iPix )
        {
            const int i = mRatio*iPix;
            QRgb color = ColorMap( A->GetImagPart(i,j), minVal, maxVal );
            painter.setPen( color );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();

    // Keep track of the extrema to allow for consistent visualization
    UpdateMinImagWindowVal( minVal );
    UpdateMaxImagWindowVal( maxVal );
}

template<typename T>
inline void 
DisplayWidget<T>::SavePng( std::string basename ) const
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWidget::SavePng");
#endif
    std::ostringstream os;
    os << basename << ".png";
    QFile file( os.str().c_str() );
    file.open( QIODevice::WriteOnly );
    pixmap_.save( &file, "PNG" );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef IO_DISPLAYWIDGET_IMPL_HPP
