/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::Syr2k
( Shape shape, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Syr2k");
    if( orientation == ConjugateTranspose )
        throw logic_error( "Syr2k accepts Normal and Transpose options." );
#endif
    if( shape == Lower && orientation == Normal )
    {
        blas::internal::Syr2kLN( alpha, A, B, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::Syr2kLT( alpha, A, B, beta, C );
    }
    else if( shape == Upper && orientation == Normal )
    {
        blas::internal::Syr2kUN( alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::Syr2kUT( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Syr2k
( Shape shape, Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Syr2k
( Shape shape, Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Syr2k
( Shape shape, Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Syr2k
( Shape shape, Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

