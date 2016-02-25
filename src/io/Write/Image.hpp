/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_WRITE_IMAGE_HPP
#define EL_WRITE_IMAGE_HPP

#ifdef EL_HAVE_QT5
# include <QFile>
# include <QImage>
# include <QPainter>
# include <QStylePainter>
#endif

namespace El {
namespace write {

#ifdef EL_HAVE_QT5
inline void
SaveQImage
( const QImage& image, string basename="matrix", 
  FileFormat format=PNG )
{
    DEBUG_ONLY(CSE cse("write::Image"))
    string filename = basename + "." + FileExtension(format);
    QFile file( filename.c_str() );
    file.open( QIODevice::WriteOnly );
    image.save( &file, QtImageFormat(format) );
}
#endif // ifdef EL_HAVE_QT5

template<typename T>
inline void
RealPartImage
( const Matrix<T>& A, string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CSE cse("write::RealPartImage"))
#ifdef EL_HAVE_QT5
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
                minVal = Min( minVal, A.GetRealPart(i,j) );
                maxVal = Max( maxVal, A.GetRealPart(i,j) );
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
#endif // ifdef EL_HAVE_QT5
}

template<typename T>
inline void
ImagPartImage
( const Matrix<T>& A, string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CSE cse("write::ImagPartImage"))
#ifdef EL_HAVE_QT5
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
                minVal = Min( minVal, A.GetImagPart(i,j) );
                maxVal = Max( maxVal, A.GetImagPart(i,j) );
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
#endif // ifdef EL_HAVE_QT5
}

template<typename Real>
inline void
Image
( const Matrix<Real>& A, string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CSE cse("write::Image"))
    RealPartImage( A, basename, format );
}

template<typename Real>
inline void
Image
( const Matrix<Complex<Real>>& A, string basename="matrix", 
  FileFormat format=PNG )
{
    DEBUG_ONLY(CSE cse("write::Image"))
    RealPartImage( A, basename+"_real", format );
    ImagPartImage( A, basename+"_imag", format );
}

} // namespace write
} // namespace El

#endif // ifndef EL_WRITE_IMAGE_HPP
