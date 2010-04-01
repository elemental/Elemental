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
#include "ElementalLAPACKInternal.h"
using namespace Elemental;

template<typename T>
void
Elemental::LAPACK::Trinv
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<T,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Trinv");
#endif
    LAPACK::Internal::TrinvVar3( shape, diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::LAPACK::Internal::TrinvVar3
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<T,MC,MR>& A  )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::TrinvVar3");
#endif
    if( shape == Lower )
        LAPACK::Internal::TrinvLVar3( diagonal, A );
    else
        LAPACK::Internal::TrinvUVar3( diagonal, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Trinv
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Internal::TrinvVar3
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<float,MC,MR>& A );

template void Elemental::LAPACK::Trinv
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

template void Elemental::LAPACK::Internal::TrinvVar3
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void Elemental::LAPACK::Trinv
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::TrinvVar3
( const Shape shape, 
  const Diagonal diagonal,
  DistMatrix<scomplex,MC,MR>& A );

template void Elemental::LAPACK::Trinv
( const Shape shape, 
  const Diagonal diagonal, 
  DistMatrix<dcomplex,MC,MR>& A );

template void Elemental::LAPACK::Internal::TrinvVar3
( const Shape shape, 
  const Diagonal diagonal,
  DistMatrix<dcomplex,MC,MR>& A );
#endif

