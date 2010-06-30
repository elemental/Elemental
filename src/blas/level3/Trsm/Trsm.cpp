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
elemental::blas::Trsm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("blas::Trsm");
#endif
    if( side == Left && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrsmLLN( diagonal, alpha, A, X );
        else
            blas::internal::TrsmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Left && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrsmLUN( diagonal, alpha, A, X );
        else
            blas::internal::TrsmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrsmRLN( diagonal, alpha, A, X );
        else
            blas::internal::TrsmRLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Upper )
    {
        if( orientation == Normal )
            blas::internal::TrsmRUN( diagonal, alpha, A, X );
        else
            blas::internal::TrsmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Trsm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::Trsm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Trsm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::Trsm
( Side side, 
  Shape shape, 
  Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

