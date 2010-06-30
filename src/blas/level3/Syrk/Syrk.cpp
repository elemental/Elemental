/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::Syrk
( Shape shape, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Syrk");
    if( orientation == ConjugateTranspose )
        throw logic_error( "Syrk accepts Normal and Transpose options." );
#endif
    if( shape == Lower && orientation == Normal )
    {
        blas::internal::SyrkLN( alpha, A, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::SyrkLT( alpha, A, beta, C );
    }
    else if( shape == Upper && orientation == Normal )
    {
        blas::internal::SyrkUN( alpha, A, beta, C );
    }
    else
    {
        blas::internal::SyrkUT( alpha, A, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Syrk
( Shape shape, Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Syrk
( Shape shape, Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Syrk
( Shape shape, Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Syrk
( Shape shape, Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

