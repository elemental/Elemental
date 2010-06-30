/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
using namespace elemental;

template<typename T>
void
elemental::blas::Symm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Symm");
#endif
    if( side == Left && shape == Lower )
    {
        blas::internal::SymmLL( alpha, A, B, beta, C );
    }
    else if( side == Left )
    {
        blas::internal::SymmLU( alpha, A, B, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::SymmRL( alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::SymmRU( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Symm
( Side side, Shape shape,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Symm
( Side side, Shape shape,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Symm
( Side side, Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Symm
( Side side, Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif
