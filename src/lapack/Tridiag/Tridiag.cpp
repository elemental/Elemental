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
elemental::lapack::Tridiag
( Shape shape, DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    if( shape == Lower )
        lapack::internal::TridiagL( A );
    else
        lapack::internal::TridiagU( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<complex<R>,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    if( shape == Lower )
        lapack::internal::TridiagL( A, t );
    else
        lapack::internal::TridiagU( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<scomplex,MD,Star>& t );

template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<dcomplex,MD,Star>& t );
#endif

