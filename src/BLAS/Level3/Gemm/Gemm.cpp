/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalBLAS_Internal.h"
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemm");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNT( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTN( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTT( orientationOfA, orientationOfB,
                                alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_A");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_A( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNT_A( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTN_A( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTT_A( orientationOfA, orientationOfB,
                                  alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_B");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_B( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNT_B( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTN_B( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTT_B( orientationOfA, orientationOfB, 
                                  alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& B,
  const T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Gemm_C");
#endif
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        BLAS::Internal::GemmNN_C( alpha, A, B, beta, C );
    }
    else if( orientationOfA == Normal )
    {
        BLAS::Internal::GemmNT_C( orientationOfB, alpha, A, B, beta, C );
    }
    else if( orientationOfB == Normal )
    {
        BLAS::Internal::GemmTN_C( orientationOfA, alpha, A, B, beta, C );
    }
    else
    {
        BLAS::Internal::GemmTT_C( orientationOfA, orientationOfB,
                                  alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& B,
  const float beta,        DistMatrix<float,MC,MR>& C );

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& B,
  const double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& B,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void Elemental::BLAS::Gemm
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_A
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_B
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );

template void Elemental::BLAS::Internal::Gemm_C
( const Orientation orientationOfA, 
  const Orientation orientationOfB,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& B,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
