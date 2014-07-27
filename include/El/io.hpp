/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IO_HPP
#define EL_IO_HPP

namespace El {

// TODO: Distributed file formats?
namespace FileFormatNS {
enum FileFormat
{
    AUTO, // Automatically detect from file extension
    ASCII,
    ASCII_MATLAB,
    BINARY,
    BINARY_FLAT,
    BMP,
    JPG,
    JPEG,
    MATRIX_MARKET,
    PNG,
    PPM,
    XBM,
    XPM,
    FileFormat_MAX // For detecting number of entries in enum
};
}
using namespace FileFormatNS;

const char* QtImageFormat( FileFormat format );
std::string FileExtension( FileFormat format );
FileFormat FormatFromExtension( const std::string ext );
FileFormat DetectFormat( const std::string filename );

std::ifstream::pos_type FileSize( std::ifstream& file );

// TODO: Many more color maps
namespace ColorMapNS {
enum ColorMap
{
    GRAYSCALE,
    GRAYSCALE_DISCRETE,
    RED_BLACK_GREEN,
    BLUE_RED
};
}
using namespace ColorMapNS;

#ifdef EL_HAVE_QT5
// Return true if Qt5 was detected, but no GUI is allowed.
// This is useful if images are to be written to file using Qt5's wrappers.
bool GuiDisabled();

// When Elemental is finalized, if no window was opened, then it must call 
// app.exit() instead
void OpenedWindow();

// For keeping track of the extreme values visualized so far
double MinRealWindowVal();
double MaxRealWindowVal();
double MinImagWindowVal();
double MaxImagWindowVal();
void UpdateMinRealWindowVal( double minVal );
void UpdateMaxRealWindowVal( double maxVal );
void UpdateMinImagWindowVal( double minVal );
void UpdateMaxImagWindowVal( double maxVal );
#endif

} // namespace El

#ifdef EL_HAVE_QT5
# include "./DisplayWidget.hpp"
# include "./DisplayWindow-premoc.hpp"
# include "./ComplexDisplayWindow-premoc.hpp"
#endif // ifdef EL_HAVE_QT5

namespace El {

// Color maps
// ==========
void SetColorMap( ColorMap colorMap );
ColorMap GetColorMap();
void SetNumDiscreteColors( Int numColors );
Int NumDiscreteColors();
#ifdef EL_HAVE_QT5
QRgb SampleColorMap( double value, double minVal, double maxVal );
#endif 

// Display
// =======
void ProcessEvents( int numMsecs );

template<typename T>
void Display( const Matrix<T>& A, std::string title="Default" );
template<typename T>
void Display( const Matrix<Complex<T>>& A, std::string title="Default" );

template<typename T,Dist U,Dist V>
void Display( const DistMatrix<T,U,V>& A, std::string title="Default" );
template<typename T,Dist U,Dist V>
void Display( const BlockDistMatrix<T,U,V>& A, std::string title="Default" );

template<typename T>
void Display( const AbstractDistMatrix<T>& AAbs, std::string title="" );
template<typename T>
void Display( const AbstractBlockDistMatrix<T>& AAbs, std::string title="" );

// Print
// =====
template<typename T>
void Print
( const Matrix<T>& A, std::string title="", std::ostream& os=std::cout );

template<typename T,Dist U,Dist V>
void Print
( const DistMatrix<T,U,V>& A,
  std::string title="", std::ostream& os=std::cout );
template<typename T,Dist U,Dist V>
void Print
( const BlockDistMatrix<T,U,V>& A,
  std::string title="", std::ostream& os=std::cout );

template<typename T>
void Print
( const AbstractDistMatrix<T>& AAbs, std::string title="",
  std::ostream& os=std::cout );
template<typename T>
void Print
( const AbstractBlockDistMatrix<T>& AAbs, std::string title="",
  std::ostream& os=std::cout );

// Utilities
// ---------
template<typename T>
void Print
( const std::vector<T>& x, std::string title="", std::ostream& os=std::cout );

// Read
// ====
template<typename T>
void Read( Matrix<T>& A, const std::string filename, FileFormat format=AUTO );

template<typename T,Dist U,Dist V>
void Read
( DistMatrix<T,U,V>& A, 
  const std::string filename, FileFormat format=AUTO, bool sequential=false );
template<typename T,Dist U,Dist V>
void Read
( BlockDistMatrix<T,U,V>& A, 
  const std::string filename, FileFormat format=AUTO, bool sequential=false );

template<typename T>
void Read
( AbstractDistMatrix<T>& A, 
  const std::string filename, FileFormat format=AUTO, bool sequential=false );
template<typename T>
void Read
( AbstractBlockDistMatrix<T>& A, 
  const std::string filename, FileFormat format=AUTO, bool sequential=false );

// Spy
// ===
template<typename T>
void Spy( const Matrix<T>& A, std::string title="Default", Base<T> tol=0 );

template<typename T,Distribution U,Distribution V>
void Spy
( const DistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 );
template<typename T,Dist U,Dist V>
void Spy
( const BlockDistMatrix<T,U,V>& A, std::string title="Default", Base<T> tol=0 );

template<typename T>
void Spy
( const AbstractDistMatrix<T>& A, std::string title="Default", Base<T> tol=0 );
template<typename T>
void Spy
( const AbstractBlockDistMatrix<T>& A,
  std::string title="Default", Base<T> tol=0 );

// Write
// =====
template<typename T>
void Write
( const Matrix<T>& A, std::string basename="matrix", FileFormat format=BINARY,
  std::string title="" );

template<typename T,Dist U,Dist V>
void Write
( const DistMatrix<T,U,V>& A, std::string basename="matrix",
  FileFormat format=BINARY, std::string title="" );
template<typename T,Dist U,Dist V>
void Write
( const BlockDistMatrix<T,U,V>& A, std::string basename="matrix",
  FileFormat format=BINARY, std::string title="" );

template<typename T>
void Write
( const AbstractDistMatrix<T>& A, std::string basename="matrix",
  FileFormat format=BINARY, std::string title="" );
template<typename T>
void Write
( const AbstractBlockDistMatrix<T>& A, std::string basename="matrix",
  FileFormat format=BINARY, std::string title="" );

} // namespace El

#endif // ifndef EL_IO_HPP
