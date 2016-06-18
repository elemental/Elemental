/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Write/Ascii.hpp"
#include "./Write/AsciiMatlab.hpp"
#include "./Write/Binary.hpp"
#include "./Write/BinaryFlat.hpp"
#include "./Write/Image.hpp"
#include "./Write/MatrixMarket.hpp"

namespace El {

template<typename T>
void Write
( const Matrix<T>& A, 
  string basename, FileFormat format, string title )
{
    DEBUG_CSE
    switch( format )
    {
    case ASCII:         write::Ascii( A, basename, title );       break;
    case ASCII_MATLAB:  write::AsciiMatlab( A, basename, title ); break;
    case BINARY:        write::Binary( A, basename );             break;
    case BINARY_FLAT:   write::BinaryFlat( A, basename );         break;
    case MATRIX_MARKET: write::MatrixMarket( A, basename );       break;
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

template<typename T>
void Write
( const AbstractDistMatrix<T>& A, 
  string basename, FileFormat format, string title )
{
    DEBUG_CSE
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Write( A.LockedMatrix(), basename, format, title );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Write( A_CIRC_CIRC.LockedMatrix(), basename, format, title );
    }
}

#define PROTO(T) \
  template void Write \
  ( const Matrix<T>& A, \
    string basename, FileFormat format, string title ); \
  template void Write \
  ( const AbstractDistMatrix<T>& A, \
    string basename, FileFormat format, string title );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
