/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_READ_HPP
#define ELEM_READ_HPP

#include "./Read/Ascii.hpp"
#include "./Read/AsciiMatlab.hpp"
#include "./Read/Binary.hpp"
#include "./Read/BinaryFlat.hpp"
#include "./Read/MatrixMarket.hpp"

namespace elem {

template<typename T>
inline void
Read( Matrix<T>& A, const std::string filename, FileFormat format=AUTO )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
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
        LogicError("Format unsupported for reading");
    }
}

template<typename T,Dist U,Dist V>
inline void
Read
( DistMatrix<T,U,V>& A, const std::string filename, FileFormat format=AUTO,
  bool sequential=false )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename ); 

    if( U == A.UGath && V == A.VGath )
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
        A = A_CIRC_CIRC;
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
            LogicError("Unsupported distributed read format"); 
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
Read
( BlockDistMatrix<T,U,V>& A, const std::string filename, FileFormat format=AUTO,
  bool sequential=false )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename ); 

    if( U == A.UGath && V == A.VGath )
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
        BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
        if( format == BINARY_FLAT )
            A_CIRC_CIRC.Resize( A.Height(), A.Width() );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
        {
            Read( A_CIRC_CIRC.Matrix(), filename, format );
            A_CIRC_CIRC.Resize
            ( A_CIRC_CIRC.Matrix().Height(), A_CIRC_CIRC.Matrix().Width() );
        }
        A_CIRC_CIRC.MakeSizeConsistent();
        A = A_CIRC_CIRC;
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
            LogicError("Unsupported distributed read format"); 
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_READ_HPP
