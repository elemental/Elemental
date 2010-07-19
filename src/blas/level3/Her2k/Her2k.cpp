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
elemental::blas::Her2k
( Shape shape, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Her2k");
    if( orientation == Transpose )
        throw logic_error
        ( "Her2k accepts Normal and ConjugateTranspose options." );
#endif
    if( shape == Lower && orientation == Normal )
    {
        blas::internal::Her2kLN( alpha, A, B, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::Her2kLC( alpha, A, B, beta, C );
    }
    else if( shape == Upper && orientation == Normal )
    {
        blas::internal::Her2kUN( alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::Her2kUC( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Her2k
( Shape shape, Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Her2k
( Shape shape, Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Her2k
( Shape shape, Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Her2k
( Shape shape, Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

