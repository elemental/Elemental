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
elemental::blas::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& x   )
{
#ifndef RELEASE
    PushCallStack("blas::Trsv");
#endif
    if( shape == Lower )
    {
        if( orientation == Normal )
            blas::internal::TrsvLN( diagonal, A, x );
        else
            blas::internal::TrsvLT( orientation, diagonal, A, x );
    }
    else
    {
        if( orientation == Normal )
            blas::internal::TrsvUN( diagonal, A, x );
        else
            blas::internal::TrsvUT( orientation, diagonal, A, x );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::blas::Trsv
( Shape shape, 
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& x );

template void
elemental::blas::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void
elemental::blas::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& x );

template void
elemental::blas::Trsv
( Shape shape,
  Orientation orientation,
  Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

