/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_GEMM_HPP
#define BLAS_GEMM_HPP

#include "./Gemm/NN.hpp"
#include "./Gemm/NT.hpp"
#include "./Gemm/TN.hpp"
#include "./Gemm/TT.hpp"

namespace elem {

namespace internal {

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
    PushCallStack("internal::LocalGemm");
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
      alpha, A.LockedLocalMatrix(), B.LockedLocalMatrix(),
      beta, C.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Gemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
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
#ifndef RELEASE
    PopCallStack();
#endif
}

namespace internal {

template<typename T>
inline void
GemmA
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmA");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNA( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTA( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNA( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTA( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmB
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmB");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNB( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTB( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNB( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTB( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmC
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmC");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        GemmNNC( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        GemmNTC( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNC( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTC( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmDot
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmDot");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
        GemmNNDot( alpha, A, B, beta, C );
    else
        throw std::logic_error("GemmDot only implemented for NN case");
    // This code will be enabled when the routines are implemented
    /*
    else if( orientationOfA == NORMAL )
    {
        GemmNTDot( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        GemmTNDot( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        GemmTTDot( orientationOfA, orientationOfB,
                                   alpha, A, B, beta, C );
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Gemm
( Orientation orientationOfA, 
  Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Gemm");
#endif
    if( orientationOfA == NORMAL && orientationOfB == NORMAL )
    {
        internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == NORMAL )
    {
        internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == NORMAL )
    {
        internal::GemmTN( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        internal::GemmTT
        ( orientationOfA, orientationOfB, alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_GEMM_HPP
