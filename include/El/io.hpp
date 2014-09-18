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
# include "El/io/DisplayWidget.hpp"
# include "El/io/DisplayWindow-premoc.hpp"
# include "El/io/ComplexDisplayWindow-premoc.hpp"
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

// Dense
// -----
template<typename Real>
void Display( const Matrix<Real>& A, std::string title="Matrix" );
template<typename Real>
void Display( const Matrix<Complex<Real>>& A, std::string title="Matrix" );
template<typename T>
void Display
( const AbstractDistMatrix<T>& AAbs, std::string title="DistMatrix" );
template<typename T>
void Display
( const AbstractBlockDistMatrix<T>& AAbs, std::string title="BlockDistMatrix" );

// Graphs and sparse matrices
// --------------------------
void Display( const Graph& graph, std::string title="Graph" );
void Display( const DistGraph& graph, std::string title="DistGraph" );

template<typename Real>
void Display
( const SparseMatrix<Real>& A, std::string title="SparseMatrix" );
template<typename Real>
void Display
( const SparseMatrix<Complex<Real>>& A, std::string title="SparseMatrix" );
template<typename Real>
void Display
( const DistSparseMatrix<Real>& A, 
  std::string title="DistSparseMatrix" );
template<typename Real>
void Display
( const DistSparseMatrix<Complex<Real>>& A, 
  std::string title="DistSparseMatrix" );

// Sparse-direct data structures
// -----------------------------
struct DistSymmInfo; // forward declaration
void DisplayLocal
( const DistSymmInfo& info, bool beforeFact, std::string title="" );

// Print
// =====

// Dense
// -----
template<typename T>
void Print
( const Matrix<T>& A, std::string title="Matrix", std::ostream& os=std::cout );
template<typename T>
void Print
( const AbstractDistMatrix<T>& AAbs, std::string title="DistMatrix",
  std::ostream& os=std::cout );
template<typename T>
void Print
( const AbstractBlockDistMatrix<T>& AAbs, std::string title="BlockDistMatrix",
  std::ostream& os=std::cout );

// Graphs and sparse matrices
// --------------------------
void Print
( const Graph& graph, std::string title="Graph", std::ostream& os=std::cout );
void Print
( const DistGraph& graph, std::string title="DistGraph", 
  std::ostream& os=std::cout );

template<typename T>
void Print
( const SparseMatrix<T>& A, std::string title="SparseMatrix", 
  std::ostream& os=std::cout );
template<typename T>
void Print
( const DistSparseMatrix<T>& A, std::string title="DistSparseMatrix",
  std::ostream& os=std::cout );

// Sparse-direct
// -------------
void PrintLocal
( const DistSymmInfo& info,
  std::string title="Local DistSymmInfo", std::ostream& os=std::cout );

// Utilities
// ---------
template<typename T>
void Print
( const std::vector<T>& x, std::string title="std::vector", 
  std::ostream& os=std::cout );

// Read
// ====
template<typename T>
void Read( Matrix<T>& A, const std::string filename, FileFormat format=AUTO );
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
void Spy( const Matrix<T>& A, std::string title="Matrix", Base<T> tol=0 );
template<typename T>
void Spy
( const AbstractDistMatrix<T>& A, std::string title="DistMatrix", 
  Base<T> tol=0 );
template<typename T>
void Spy
( const AbstractBlockDistMatrix<T>& A,
  std::string title="BlockDistMatrix", Base<T> tol=0 );

// Write
// =====
template<typename T>
void Write
( const Matrix<T>& A, std::string basename="Matrix", FileFormat format=BINARY,
  std::string title="" );
template<typename T>
void Write
( const AbstractDistMatrix<T>& A, std::string basename="DistMatrix",
  FileFormat format=BINARY, std::string title="" );
template<typename T>
void Write
( const AbstractBlockDistMatrix<T>& A, std::string basename="BlockDistMatrix",
  FileFormat format=BINARY, std::string title="" );

} // namespace El

#endif // ifndef EL_IO_HPP
