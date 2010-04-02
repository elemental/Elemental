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
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Trsv
( const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& x   )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trsv");
#endif
    if( shape == Lower )
    {
        if( orientation == Normal )
            BLAS::Internal::TrsvLN( diagonal, A, x );
        else
            BLAS::Internal::TrsvLT( orientation, diagonal, A, x );
    }
    else
    {
        if( orientation == Normal )
            BLAS::Internal::TrsvUN( diagonal, A, x );
        else
            BLAS::Internal::TrsvUT( orientation, diagonal, A, x );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
Elemental::BLAS::Trsv
( const Shape shape, 
  const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& x );

template void
Elemental::BLAS::Trsv
( const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void
Elemental::BLAS::Trsv
( const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& x );

template void
Elemental::BLAS::Trsv
( const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

