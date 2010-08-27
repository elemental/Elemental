/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
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
                lapack::internal::UTLLC( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTLUN( offset, H, A );
            else
                lapack::internal::UTLUC( offset, H, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTRLN( offset, H, A );
            else
                lapack::internal::UTRLC( offset, H, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTRUN( offset, H, A );
            else
                lapack::internal::UTRUC( offset, H, A );
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
                lapack::internal::UTLLC( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTLUN( offset, H, t, A );
            else
                lapack::internal::UTLUC( offset, H, t, A );
        }
    }
    else
    {
        if( shape == Lower )
        {
            if( orientation == Normal )
                lapack::internal::UTRLN( offset, H, t, A );
            else
                lapack::internal::UTRLC( offset, H, t, A );
        }
        else
        {
            if( orientation == Normal )
                lapack::internal::UTRUN( offset, H, t, A );
            else
                lapack::internal::UTRUC( offset, H, t, A );
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

