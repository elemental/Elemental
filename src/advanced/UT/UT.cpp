/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/advanced_internal.hpp"
using namespace elemental;
using namespace std;

template<typename R> // representation of a real number
void
elemental::advanced::UT
( Side side,
  Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::UT");
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
                advanced::internal::UTLLN( offset, H, A );
            else
                advanced::internal::UTLLC( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                advanced::internal::UTLUN( offset, H, A );
            else
                advanced::internal::UTLUC( offset, H, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                advanced::internal::UTRLN( offset, H, A );
            else
                advanced::internal::UTRLC( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                advanced::internal::UTRUN( offset, H, A );
            else
                advanced::internal::UTRUC( offset, H, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::advanced::UT
( Side side,
  Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<complex<R>,MC,MR  >& H, 
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::UT");
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
                advanced::internal::UTLLN( offset, H, t, A );
            else
                advanced::internal::UTLLC( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                advanced::internal::UTLUN( offset, H, t, A );
            else
                advanced::internal::UTLUC( offset, H, t, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                advanced::internal::UTRLN( offset, H, t, A );
            else
                advanced::internal::UTRLC( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                advanced::internal::UTRUN( offset, H, t, A );
            else
                advanced::internal::UTRUC( offset, H, t, A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> // representation of a real number
void
elemental::advanced::UT
( Side side,
  Shape shape, 
  Orientation orientation,
  int offset,
  const DistMatrix<complex<R>,MC,  MR  >& H, 
  const DistMatrix<complex<R>,Star,Star>& t,
        DistMatrix<complex<R>,MC,  MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::UT");
    if( orientation == Transpose )
        throw logic_error
              ( "Only Normal and ConjugateTranspose UT transform applications "
                "are written." );
#endif
    DistMatrix<complex<R>,MD,Star> tDiag(A.Grid());
    tDiag.AlignWithDiag( A, offset );
    tDiag = t;
    advanced::UT( side, shape, orientation, offset, H, tDiag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif // WITHOUT_COMPLEX

template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<scomplex,MC,  MR  >& H,
  const DistMatrix<scomplex,MD,  Star>& t,
        DistMatrix<scomplex,MC,  MR  >& A );
template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<scomplex,MC,  MR  >& H,
  const DistMatrix<scomplex,Star,Star>& t,
        DistMatrix<scomplex,MC,  MR  >& A );

template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<dcomplex,MC,  MR  >& H,
  const DistMatrix<dcomplex,MD,  Star>& t,
        DistMatrix<dcomplex,MC,  MR  >& A );
template void elemental::advanced::UT
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  const DistMatrix<dcomplex,MC,  MR  >& H,
  const DistMatrix<dcomplex,Star,Star>& t,
        DistMatrix<dcomplex,MC,  MR  >& A );
#endif

