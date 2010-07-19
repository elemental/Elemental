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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::Herk
( Shape shape, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Herk");
    if( A.GetGrid() != C.GetGrid() )
        throw logic_error( "A and C must be distributed over the same grid." );
    if( orientation == Transpose )
        throw logic_error
        ( "Herk accepts Normal and ConjugateTranspose options." );
#endif
    if( shape == Lower && orientation == Normal )
    {
        blas::internal::HerkLN( alpha, A, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::HerkLC( alpha, A, beta, C );
    }
    else if( shape == Upper && orientation == Normal )
    {
        blas::internal::HerkUN( alpha, A, beta, C );
    }
    else
    {
        blas::internal::HerkUC( alpha, A, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Herk
( Shape shape, Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Herk
( Shape shape, Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Herk
( Shape shape, Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Herk
( Shape shape, Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

