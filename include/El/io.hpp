/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IO_HPP
#define EL_IO_HPP

namespace El {

const char* QtImageFormat( FileFormat format );
string FileExtension( FileFormat format );
FileFormat FormatFromExtension( const string ext );
FileFormat DetectFormat( const string filename );

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
void Display( const Matrix<Real>& A, string title="Matrix" );
template<typename Real>
void Display( const Matrix<Complex<Real>>& A, string title="Matrix" );
template<typename T>
void Display( const AbstractDistMatrix<T>& A, string title="DistMatrix" );
template<typename T>
void Display( const DistMultiVec<T>& X, string title="DistMultiVec" );

// Graphs and sparse matrices
// --------------------------
void Display( const Graph& graph, string title="Graph" );
void Display( const DistGraph& graph, string title="DistGraph" );

template<typename Real>
void Display
( const SparseMatrix<Real>& A, string title="SparseMatrix" );
template<typename Real>
void Display
( const SparseMatrix<Complex<Real>>& A, string title="SparseMatrix" );
template<typename T>
void Display( const DistSparseMatrix<T>& A, string title="DistSparseMatrix" );

// Sparse-direct data structures
// -----------------------------
namespace ldl { struct DistNodeInfo; }
void DisplayLocal
( const ldl::DistNodeInfo& info, bool beforeFact, string title="" );

// Print
// =====

// Dense
// -----
template<typename T>
void Print( const Matrix<T>& A, string title="Matrix", ostream& os=cout );
template<typename T>
void Print
( const AbstractDistMatrix<T>& A, string title="DistMatrix", ostream& os=cout );
template<typename T>
void Print
( const DistMultiVec<T>& X, string title="DistMultiVec", ostream& os=cout );

// Graphs and sparse matrices
// --------------------------
void Print( const Graph& graph, string title="Graph", ostream& os=cout );
void Print
( const DistGraph& graph, string title="DistGraph", ostream& os=cout );

template<typename T>
void Print
( const SparseMatrix<T>& A, string title="SparseMatrix", ostream& os=cout );
template<typename T>
void Print
( const DistSparseMatrix<T>& A, string title="DistSparseMatrix",
  ostream& os=cout );

// Sparse-direct
// -------------
void PrintLocal
( const ldl::DistNodeInfo& info,
  string title="Local ldl::DistNodeInfo", ostream& os=cout );

// Utilities
// ---------
template<typename T>
void Print( const vector<T>& x, string title="vector", ostream& os=cout );

// Read
// ====
template<typename T>
void Read( Matrix<T>& A, const string filename, FileFormat format=AUTO );
template<typename T>
void Read
( AbstractDistMatrix<T>& A, 
  const string filename, FileFormat format=AUTO, bool sequential=false );

// Spy
// ===
template<typename T>
void Spy( const Matrix<T>& A, string title="Matrix", Base<T> tol=0 );
template<typename T>
void Spy
( const AbstractDistMatrix<T>& A, string title="DistMatrix", Base<T> tol=0 );

// Write
// =====
template<typename T>
void Write
( const Matrix<T>& A, string basename="Matrix", FileFormat format=BINARY,
  string title="" );
template<typename T>
void Write
( const AbstractDistMatrix<T>& A, string basename="DistMatrix",
  FileFormat format=BINARY, string title="" );

} // namespace El

#endif // ifndef EL_IO_HPP
