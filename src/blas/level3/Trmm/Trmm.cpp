/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace elemental;

template<typename T>
void
elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("blas::Trmm");
#endif
    if( side == Left && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrmmLLN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Left && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrmmLUN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrmmRLN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmRLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrmmRUN( diagonal, alpha, A, X );
        else
            blas::internal::TrmmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::Trmm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

