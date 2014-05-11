/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_WRITE_IMAGE_HPP
#define ELEM_WRITE_IMAGE_HPP

#ifdef ELEM_HAVE_QT5
# include <QFile>
# include <QImage>
# include <QPainter>
# include <QStylePainter>
#endif

namespace elem {
namespace write {

#ifdef ELEM_HAVE_QT5
inline void
SaveQImage
( const QImage& image, std::string basename="matrix", 
  FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::Image"))
    std::string filename = basename + "." + FileExtension(format);
    QFile file( filename.c_str() );
    file.open( QIODevice::WriteOnly );
    image.save( &file, QtImageFormat(format) );
}
#endif // ifdef ELEM_HAVE_QT5

template<typename T>
inline void
RealPartImage
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::RealPartImage"))
#ifdef ELEM_HAVE_QT5
    typedef Base<T> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    // Compute the maximum and minimum values
    Real minVal=0, maxVal=0; 
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A.GetRealPart( 0, 0 );
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<m; ++i )
            {
                minVal = std::min( minVal, A.GetRealPart(i,j) );
                maxVal = std::max( maxVal, A.GetRealPart(i,j) );
            }
        }
    }

    // TODO: Parameterize these instead
    const Int mPix = 2*m;
    const Int nPix = 2*n;
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    QImage image( nPix, mPix, QImage::Format_RGB32 );
    for( Int jPix=0; jPix<nPix; ++jPix )
    {
        const Int j = nRatio*jPix;
        for( Int iPix=0; iPix<mPix; ++iPix )
        {
            const Int i = mRatio*iPix;
            QRgb color = SampleColorMap( A.GetRealPart(i,j), minVal, maxVal );
            image.setPixel( jPix, iPix, color );
        }
    }

    SaveQImage( image, basename, format );
#else
    LogicError("Qt5 not available");
#endif // ifdef ELEM_HAVE_QT5
}

template<typename T>
inline void
ImagPartImage
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::ImagPartImage"))
#ifdef ELEM_HAVE_QT5
    typedef Base<T> Real;
    const Int m = A.Height();
    const Int n = A.Width();

    // Compute the maximum and minimum values
    Real minVal=0, maxVal=0; 
    if( m != 0 && n != 0 )
    {
        minVal = maxVal = A.GetImagPart( 0, 0 );
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<m; ++i )
            {
                minVal = std::min( minVal, A.GetImagPart(i,j) );
                maxVal = std::max( maxVal, A.GetImagPart(i,j) );
            }
        }
    }

    // TODO: Parameterize these instead
    const Int mPix = 2*m;
    const Int nPix = 2*n;
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    QImage image( nPix, mPix, QImage::Format_RGB32 );
    for( Int jPix=0; jPix<nPix; ++jPix )
    {
        const Int j = nRatio*jPix;
        for( Int iPix=0; iPix<mPix; ++iPix )
        {
            const Int i = mRatio*iPix;
            QRgb color = SampleColorMap( A.GetImagPart(i,j), minVal, maxVal );
            image.setPixel( jPix, iPix, color );
        }
    }

    SaveQImage( image, basename, format );
#else
    LogicError("Qt5 not available");
#endif // ifdef ELEM_HAVE_QT5
}

template<typename Real>
inline void
Image
( const Matrix<Real>& A, std::string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::Image"))
    RealPartImage( A, basename, format );
}

template<typename Real>
inline void
Image
( const Matrix<Complex<Real>>& A, std::string basename="matrix", 
  FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::Image"))
    RealPartImage( A, basename+"_real", format );
    ImagPartImage( A, basename+"_imag", format );
}

} // namespace write
} // namespace elem

#endif // ifndef ELEM_WRITE_IMAGE_HPP
