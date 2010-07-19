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
elemental::blas::Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A, 
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Gemv");
#endif
    if( orientation == Normal )
        blas::internal::GemvN( alpha, A, x, beta, y );
    else
        blas::internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Gemv
( Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& x,
  float beta,        DistMatrix<float,MC,MR>& y );

template void elemental::blas::Gemv
( Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& x,
  double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Gemv
( Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& x,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void elemental::blas::Gemv
( Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& x,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif

