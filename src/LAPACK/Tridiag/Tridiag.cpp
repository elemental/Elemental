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
Elemental::LAPACK::Tridiag
( const Shape shape,
  DistMatrix<T,MC,MR  >& A,
  DistMatrix<T,MD,Star>& d,
  DistMatrix<T,MD,Star>& e,
  DistMatrix<T,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Tridiag");
#endif
    if( shape == Lower )
        LAPACK::Internal::TridiagL( A, d, e, t );
    else
        LAPACK::Internal::TridiagU( A, d, e, t );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Tridiag
( const Shape shape, 
  DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MD,Star>& d,
  DistMatrix<float,MD,Star>& e,
  DistMatrix<float,MD,Star>& t );

template void Elemental::LAPACK::Tridiag
( const Shape shape, 
  DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MD,Star>& d,
  DistMatrix<double,MD,Star>& e,
  DistMatrix<double,MD,Star>& t );

