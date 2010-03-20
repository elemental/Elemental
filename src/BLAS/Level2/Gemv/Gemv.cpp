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
#include "ElementalBLAS_Internal.h"
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Gemv
( const Orientation orientation,
  const T alpha, const DistMatrix<T,MC,MR>& A,
                 const DistMatrix<T,MC,MR>& x,
  const T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemv");
#endif
    if( orientation == Normal )
        BLAS::Internal::GemvN( alpha, A, x, beta, y );
    else
        BLAS::Internal::GemvT( orientation, alpha, A, x, beta, y );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const float alpha, const DistMatrix<float,MC,MR>& A,
                     const DistMatrix<float,MC,MR>& x,
  const float beta,        DistMatrix<float,MC,MR>& y );

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const double alpha, const DistMatrix<double,MC,MR>& A,
                      const DistMatrix<double,MC,MR>& x,
  const double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                        const DistMatrix<scomplex,MC,MR>& x,
  const scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void Elemental::BLAS::Gemv
( const Orientation orientation,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                        const DistMatrix<dcomplex,MC,MR>& x,
  const dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif
