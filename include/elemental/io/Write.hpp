/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_WRITE_HPP
#define ELEM_IO_WRITE_HPP

#include "elemental/io/Print.hpp"

#ifdef HAVE_QT5
# include <QFile>
# include <QImage>
# include <QPainter>
# include <QStylePainter>
#endif

namespace elem {
namespace write {

template<typename T>
inline void
Ascii( const Matrix<T>& A, std::string basename="matrix", std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("write::Ascii"))
    std::string filename = basename + "." + FileExtension(ASCII);
    std::ofstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    file.setf( std::ios::scientific );
    Print( A, title, file );
}

template<typename T>
inline void
AsciiMatlab
( const Matrix<T>& A, std::string basename="matrix", 
  std::string title="matrix" )
{
    DEBUG_ONLY(CallStackEntry cse("write::AsciiMatlab"))
    // Empty titles are not legal
    if( title == "" )
        title = "matrix";

    std::string filename = basename + "." + FileExtension(ASCII_MATLAB);
    std::ofstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    file.setf( std::ios::scientific );
    file << title << " = [\n";
    Print( A, "", file );
    file << "];\n";
}

template<typename T>
inline void
BinaryFlat( const Matrix<T>& A, std::string basename="matrix" )
{
    DEBUG_ONLY(CallStackEntry cse("write::BinaryFlat"))
    
    std::string filename = basename + "." + FileExtension(BINARY_FLAT);
    std::ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    if( A.Height() == A.LDim() )
        file.write( (char*)A.LockedBuffer(), A.Height()*A.Width()*sizeof(T) );
    else
        for( Int j=0; j<A.Width(); ++j )
            file.write( (char*)A.LockedBuffer(0,j), A.Height()*sizeof(T) );
}

template<typename T>
inline void
Binary( const Matrix<T>& A, std::string basename="matrix" )
{
    DEBUG_ONLY(CallStackEntry cse("write::Binary"))
    
    std::string filename = basename + "." + FileExtension(BINARY);
    std::ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    file << A.Height() << A.Width();
    if( A.Height() == A.LDim() )
        file.write( (char*)A.LockedBuffer(), A.Height()*A.Width()*sizeof(T) );
    else
        for( Int j=0; j<A.Width(); ++j )
            file.write( (char*)A.LockedBuffer(0,j), A.Height()*sizeof(T) );
}

#ifdef HAVE_QT5
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
#endif // ifdef HAVE_QT5

template<typename T>
inline void
RealPartImage
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::RealPartImage"))
#ifdef HAVE_QT5
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
    const Int mPix = Max( 500, 2*m );
    const Int nPix = Max( 500, 2*n );
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
#endif // ifdef HAVE_QT5
}

template<typename T>
inline void
ImagPartImage
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=PNG )
{
    DEBUG_ONLY(CallStackEntry cse("write::ImagPartImage"))
#ifdef HAVE_QT5
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
    const Int mPix = Max( 500, 2*m );
    const Int nPix = Max( 500, 2*n );
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
#endif // ifdef HAVE_QT5
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

template<typename T>
inline void
Write
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=BINARY, 
  std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("Write"))
    switch( format )
    {
    case ASCII:        write::Ascii( A, basename, title );       break;
    case ASCII_MATLAB: write::AsciiMatlab( A, basename, title ); break;
    case BINARY:       write::Binary( A, basename );             break;
    case BINARY_FLAT:  write::BinaryFlat( A, basename );         break;
    case BMP:
    case JPG:
    case JPEG:
    case PNG:
    case PPM:
    case XBM:
    case XPM:
        write::Image( A, basename, format ); break;
    default:
        LogicError("Invalid file format");
    }
}

template<typename T,Distribution U,Distribution V>
inline void
Write
( const DistMatrix<T,U,V>& A, std::string basename="matrix", 
  FileFormat format=BINARY, std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("Write"))
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
    if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
        Write( A_CIRC_CIRC.LockedMatrix(), basename, format, title );
}

// If already in [* ,* ] or [o ,o ] distributions, no copy is needed

template<typename T>
inline void
Write
( const DistMatrix<T,STAR,STAR>& A, std::string basename="matrix", 
  FileFormat format=BINARY, std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("Write"))
    if( A.Grid().VCRank() == 0 )
        Write( A.LockedMatrix(), basename, format, title );
}
template<typename T>
inline void
Write
( const DistMatrix<T,CIRC,CIRC>& A, std::string basename="matrix", 
  FileFormat format=BINARY, std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("Write"))
    if( A.CrossRank() == A.Root() )
        Write( A.LockedMatrix(), basename, format, title );
}

} // namespace elem

#endif // ifndef ELEM_IO_WRITE_HPP
