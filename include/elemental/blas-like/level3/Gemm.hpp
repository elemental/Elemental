/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_GEMM_HPP
#define ELEM_BLAS_GEMM_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T,Distribution AColDist,Distribution ARowDist,
                    Distribution BColDist,Distribution BRowDist,
                    Distribution CColDist,Distribution CRowDist>
inline void LocalGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,BColDist,BRowDist>& B,
  T beta,        DistMatrix<T,CColDist,CRowDist>& C )
{
#ifndef RELEASE
    CallStackEntry entry("LocalGemm");
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( AColDist != CColDist ||
            ARowDist != BColDist ||
            BRowDist != CRowDist )
            throw std::logic_error("C[X,Y] = A[X,Z] B[Z,Y]");
        if( A.ColAlignment() != C.ColAlignment() )
            throw std::logic_error("A's cols must align with C's rows");
        if( A.RowAlignment() != B.ColAlignment() )
            throw std::logic_error("A's rows must align with B's cols");
        if( B.RowAlignment() != C.RowAlignment() )
            throw std::logic_error("B's rows must align with C's rows");
        if( A.Height() != C.Height() ||
            A.Width() != B.Height() ||
            B.Width() != C.Width() )
        {
            std::ostringstream msg;
            msg << "Nonconformal LocalGemmNN:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
                << "  C ~ " << C.Height() << " x " << C.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else if( orientationOfA == NORMAL )
    {
        if( AColDist != CColDist ||
            ARowDist != BRowDist ||
            BColDist != CRowDist )
            throw std::logic_error("C[X,Y] = A[X,Z] (B[Y,Z])^(T/H)");
        if( A.ColAlignment() != C.ColAlignment() )
            throw std::logic_error("A's cols must align with C's rows");
        if( A.RowAlignment() != B.RowAlignment() )
            throw std::logic_error("A's rows must align with B's rows");
        if( B.ColAlignment() != C.RowAlignment() )
            throw std::logic_error("B's cols must align with C's rows");
        if( A.Height() != C.Height() ||
            A.Width() != B.Width() ||
            B.Height() != C.Width() )
        {
            std::ostringstream msg;
            msg << "Nonconformal LocalGemmNT:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
                << "  C ~ " << C.Height() << " x " << C.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else if( orientationOfB == NORMAL )
    {
        if( ARowDist != CColDist ||
            AColDist != BColDist ||
            BRowDist != CRowDist )
            throw std::logic_error("C[X,Y] = (A[Z,X])^(T/H) B[Z,Y]");
        if( A.RowAlignment() != C.ColAlignment() )
            throw std::logic_error("A's rows must align with C's cols");
        if( A.ColAlignment() != B.ColAlignment() )
            throw std::logic_error("A's cols must align with B's cols");
        if( B.RowAlignment() != C.RowAlignment() )
            throw std::logic_error("B's rows must align with C's rows");
        if( A.Width() != C.Height() ||
            A.Height() != B.Height() ||
            B.Width() != C.Width() )
        {
            std::ostringstream msg;
            msg << "Nonconformal LocalGemmTN:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
                << "  C ~ " << C.Height() << " x " << C.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
    else
    {
        if( ARowDist != CColDist ||
            AColDist != BRowDist ||
            BColDist != CRowDist )
            throw std::logic_error("C[X,Y] = (A[Z,X])^(T/H) (B[Y,Z])^(T/H)");
        if( A.RowAlignment() != C.ColAlignment() )
            throw std::logic_error("A's rows must align with C's cols");
        if( A.ColAlignment() != B.RowAlignment() )
            throw std::logic_error("A's cols must align with B's rows");
        if( B.ColAlignment() != C.RowAlignment() )
            throw std::logic_error("B's cols must align with C's rows");
        if( A.Width() != C.Height() ||
            A.Height() != B.Width() ||
            B.Height() != C.Width() )
        {
            std::ostringstream msg;
            msg << "Nonconformal LocalGemmTT:\n"
                << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
                << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
                << "  C ~ " << C.Height() << " x " << C.Width();
            throw std::logic_error( msg.str().c_str() );
        }
    }
#endif
    Gemm
    ( orientationOfA , orientationOfB,
      alpha, A.LockedMatrix(), B.LockedMatrix(), beta, C.Matrix() );
}

template<typename T,Distribution AColDist,Distribution ARowDist,
                    Distribution BColDist,Distribution BRowDist,
                    Distribution CColDist,Distribution CRowDist>
inline void LocalGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,BColDist,BRowDist>& B,
                 DistMatrix<T,CColDist,CRowDist>& C )
{
#ifndef RELEASE
    CallStackEntry entry("LocalGemm");
#endif
    const int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    LocalGemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C );
}

} // namespace elem

#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace elem {

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Gemm");
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height() )
            throw std::logic_error("Nonconformal GemmNN");
    }
    else if( orientationOfA == NORMAL )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width() )
            throw std::logic_error("Nonconformal GemmN(T/C)");
    }
    else if( orientationOfB == NORMAL )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height() )
            throw std::logic_error("Nonconformal Gemm(T/C)N");
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width() )
            throw std::logic_error("Nonconformal Gemm(T/C)(T/C)");
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == NORMAL ? A.Width() : A.Height() );
    if( k != 0 )
    {
        blas::Gemm
        ( transA, transB, m, n, k,
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Scale( beta, C );
    }
}

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Gemm");
#endif
    const int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C );
}

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Gemm");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        gemm::SUMMA_NN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        gemm::SUMMA_NT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        gemm::SUMMA_TN( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        gemm::SUMMA_TT( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
}

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry entry("Gemm");
#endif
    const int m = ( orientationOfA==NORMAL ? A.Height() : A.Width() );
    const int n = ( orientationOfB==NORMAL ? B.Width() : B.Height() );
    Zeros( C, m, n );
    Gemm( orientationOfA, orientationOfB, alpha, A, B, T(0), C );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_GEMM_HPP
