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
using namespace elemental;

template<typename T>
void
elemental::blas::Hemm
(  Side side, Shape shape,
   T alpha, const DistMatrix<T,MC,MR>& A,
            const DistMatrix<T,MC,MR>& B,
   T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::Hemm");
#endif
    if( side == Left && shape == Lower )
    {
        blas::internal::HemmLL( alpha, A, B, beta, C );
    }
    else if( side == Left )
    {
        blas::internal::HemmLU( alpha, A, B, beta, C );
    }
    else if( shape == Lower )
    {
        blas::internal::HemmRL( alpha, A, B, beta, C );
    }
    else
    {
        blas::internal::HemmRU( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Hemm
( Side side, Shape shape,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::Hemm
( Side side, Shape shape,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Hemm
( Side side, Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::Hemm
( Side side, Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

