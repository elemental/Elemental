/*
   This file is part of elemental, a library for distributed-memory dense
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the
   file LICENSE.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

template<typename T>
void
elemental::lapack::UT
( Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<T,MC,MR>& H, 
        DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::UT");
    if( orientation == Transpose )
        throw logic_error
              ( "Only Normal and ConjugateTranspose UT transform applications "
                "are written." );
#endif
    if( shape == Lower )
    {
        if( orientation == Normal )
            lapack::internal::UTLN( offset, H, A );
        else
            lapack::internal::UTLH( offset, H, A );
    }
    else
    {
        if( orientation == Normal )
            lapack::internal::UTUN( offset, H, A );
        else
            lapack::internal::UTUH( offset, H, A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::UT
( Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::lapack::UT
( Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::UT
( Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<scomplex,MC,MR>& H,
        DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::UT
( Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<dcomplex,MC,MR>& H,
        DistMatrix<dcomplex,MC,MR>& A );
#endif

