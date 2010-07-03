/*
   This file is part of Elemental, a library for distributed-memory dense
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

template<typename R>
void
elemental::lapack::UT
( Side side,
  Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::UT");
    if( orientation == Transpose )
        throw logic_error
              ( "Only Normal and ConjugateTranspose UT transform applications "
                "are written." );
#endif
    if( side == Left )
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTLLN( offset, H, A );
            else
                lapack::internal::UTLLH( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTLUN( offset, H, A );
            else
                lapack::internal::UTLUH( offset, H, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTRLN( offset, H, A );
            else
                lapack::internal::UTRLH( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTRUN( offset, H, A );
            else
                lapack::internal::UTRUH( offset, H, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::UT
( Side side,
  Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<complex<R>,MC,MR  >& H, 
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("lapack::UT");
    if( orientation == Transpose )
        throw logic_error
              ( "Only Normal and ConjugateTranspose UT transform applications "
                "are written." );
#endif
    if( side == Left )
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTLLN( offset, H, t, A );
            else
                lapack::internal::UTLLH( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTLUN( offset, H, t, A );
            else
                lapack::internal::UTLUH( offset, H, t, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTRLN( offset, H, t, A );
            else
                lapack::internal::UTRLH( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTRUN( offset, H, t, A );
            else
                lapack::internal::UTRUH( offset, H, t, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::lapack::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<scomplex,MC,MR  >& H,
  const DistMatrix<scomplex,MD,Star>& t,
        DistMatrix<scomplex,MC,MR  >& A );

template void elemental::lapack::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<dcomplex,MC,MR  >& H,
  const DistMatrix<dcomplex,MD,Star>& t,
        DistMatrix<dcomplex,MC,MR  >& A );
#endif

