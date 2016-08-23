/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Read/Ascii.hpp"
#include "./Read/AsciiMatlab.hpp"
#include "./Read/Binary.hpp"
#include "./Read/BinaryFlat.hpp"
#include "./Read/MatrixMarket.hpp"

namespace El {

template<typename T>
void Read( Matrix<T>& A, const string filename, FileFormat format )
{
    DEBUG_CSE
    if( format == AUTO )
        format = DetectFormat( filename );

    switch( format )
    {
    case ASCII:
        read::Ascii( A, filename );
        break;
    case ASCII_MATLAB:
        read::AsciiMatlab( A, filename );
        break;
    case BINARY:
        read::Binary( A, filename );
        break;
    case BINARY_FLAT:
        read::BinaryFlat( A, A.Height(), A.Width(), filename );
        break;
    case MATRIX_MARKET:
        read::MatrixMarket( A, filename );
        break;
    default:
        LogicError("Format unsupported for reading a Matrix");
    }
}

template<typename T>
void Read
( AbstractDistMatrix<T>& A, const string filename, FileFormat format,
  bool sequential )
{
    DEBUG_CSE
    if( format == AUTO )
        format = DetectFormat( filename ); 

    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
        {
            Read( A.Matrix(), filename, format );
            A.Resize( A.Matrix().Height(), A.Matrix().Width() );
        }
        A.MakeSizeConsistent();
    }
    else if( sequential )
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
        if( format == BINARY_FLAT )
            A_CIRC_CIRC.Resize( A.Height(), A.Width() );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
        {
            Read( A_CIRC_CIRC.Matrix(), filename, format );
            A_CIRC_CIRC.Resize
            ( A_CIRC_CIRC.Matrix().Height(), A_CIRC_CIRC.Matrix().Width() );
        }
        A_CIRC_CIRC.MakeSizeConsistent();
        Copy( A_CIRC_CIRC, A );
    }
    else
    {
        switch( format )
        {
        case ASCII:
            read::Ascii( A, filename );
            break;
        case ASCII_MATLAB:
            read::AsciiMatlab( A, filename );
            break;
        case BINARY:
            read::Binary( A, filename );
            break;
        case BINARY_FLAT:
            read::BinaryFlat( A, A.Height(), A.Width(), filename );
            break;
        case MATRIX_MARKET:
            read::MatrixMarket( A, filename );
            break;
        default:
            LogicError("Format unsupported for reading a DistMatrix"); 
        }
    }
}

template<typename T>
void Read( SparseMatrix<T>& A, const string filename, FileFormat format )
{
    DEBUG_CSE
    if( format == AUTO )
        format = DetectFormat( filename );

    switch( format )
    {
    case MATRIX_MARKET:
        read::MatrixMarket( A, filename );
        break;
    default:
        LogicError("Format unsupported for reading a SparseMatrix");
    }
}

template<typename T>
void Read( DistSparseMatrix<T>& A, const string filename, FileFormat format )
{
    DEBUG_CSE
    if( format == AUTO )
        format = DetectFormat( filename );

    switch( format )
    {
    case MATRIX_MARKET:
        read::MatrixMarket( A, filename );
        break;
    default:
        LogicError("Format unsupported for reading a DistSparseMatrix");
    }
}

#define PROTO(T) \
  template void Read \
  ( Matrix<T>& A, const string filename, FileFormat format ); \
  template void Read \
  ( AbstractDistMatrix<T>& A, const string filename, \
    FileFormat format, bool sequential ); \
  template void Read \
  ( SparseMatrix<T>& A, const string filename, FileFormat format ); \
  template void Read \
  ( DistSparseMatrix<T>& A, const string filename, FileFormat format );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
