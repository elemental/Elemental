/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El/io.hpp"

#ifdef EL_HAVE_QT5

#include <QFile>
#include <QPainter>
#include <QPixmap>
#include <QStylePainter>

#include "El/io/DisplayWidget.hpp"

namespace El {

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
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::paintEvent"))
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap_ );
}

template<typename T>
void DisplayWidget<T>::DisplayReal( const Matrix<T>* A )
{
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::DisplayReal"))
    typedef Base<T> Real;
    const Int m = A->Height();
    const Int n = A->Width();
    if( m == 0 || n == 0 )
        return;

    // Compute the range of the real values in A
    Real minVal=0, maxVal=0;
    minVal = maxVal = A->GetRealPart( 0, 0 );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            minVal = std::min( minVal, A->GetRealPart(i,j) );
            maxVal = std::max( maxVal, A->GetRealPart(i,j) );
        }
    }

    DisplayReal( A, minVal, maxVal );
}

template<typename T>
void DisplayWidget<T>::DisplayReal
( const Matrix<T>* A, Base<T> minVal, Base<T> maxVal )
{
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::DisplayReal"))
    const Int m = A->Height();
    const Int n = A->Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Parameterize these instead
    const Int mPix = 2*m;
    const Int nPix = 2*n;
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( Int jPix=0; jPix<nPix; ++jPix ) 
    {
        const Int j = nRatio*jPix;
        for( Int iPix=0; iPix<mPix; ++iPix )
        {
            const Int i = mRatio*iPix;
            QRgb color = SampleColorMap( A->GetRealPart(i,j), minVal, maxVal );
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
void DisplayWidget<T>::DisplayImag( const Matrix<T>* A )
{
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::DisplayImag"))
    typedef Base<T> Real;
    const Int m = A->Height();
    const Int n = A->Width();
    if( m == 0 || n == 0 )
        return;

    // Compute the range of the real values in A
    Real minVal=0, maxVal=0;
    minVal = maxVal = A->GetImagPart( 0, 0 );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            minVal = std::min( minVal, A->GetImagPart(i,j) );
            maxVal = std::max( maxVal, A->GetImagPart(i,j) );
        }
    }

    DisplayImag( A, minVal, maxVal );
}

template<typename T>
void DisplayWidget<T>::DisplayImag
( const Matrix<T>* A, Base<T> minVal, Base<T> maxVal )
{
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::DisplayImag"))
    const Int m = A->Height();
    const Int n = A->Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Parameterize these instead
    const Int mPix = 2*m;
    const Int nPix = 2*n;
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( Int jPix=0; jPix<nPix; ++jPix ) 
    {
        const Int j = nRatio*jPix;
        for( Int iPix=0; iPix<mPix; ++iPix )
        {
            const Int i = mRatio*iPix;
            QRgb color = SampleColorMap( A->GetImagPart(i,j), minVal, maxVal );
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
void DisplayWidget<T>::SavePng( std::string basename ) const
{
    DEBUG_ONLY(CallStackEntry cse("DisplayWidget::SavePng"))
    std::string filename = basename + ".png";
    QFile file( filename.c_str() );
    file.open( QIODevice::WriteOnly );
    pixmap_.save( &file, "PNG" );
}

#define PROTO(T) template class DisplayWidget<T>

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El

#endif // ifdef EL_HAVE_QT5
