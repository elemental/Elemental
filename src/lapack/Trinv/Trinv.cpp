/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;

template<typename T>
void
elemental::lapack::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<T,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("lapack::Trinv");
#endif
    lapack::internal::TrinvVar3( shape, diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<T,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvVar3");
#endif
    if( shape == Lower )
        lapack::internal::TrinvLVar3( diagonal, A );
    else
        lapack::internal::TrinvUVar3( diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

template void elemental::lapack::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal,
  DistMatrix<scomplex,MC,MR>& A );

template void elemental::lapack::Trinv
( Shape shape, 
  Diagonal diagonal, 
  DistMatrix<dcomplex,MC,MR>& A );

template void elemental::lapack::internal::TrinvVar3
( Shape shape, 
  Diagonal diagonal,
  DistMatrix<dcomplex,MC,MR>& A );
#endif

