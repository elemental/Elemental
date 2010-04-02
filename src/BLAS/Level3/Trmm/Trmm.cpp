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
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const T alpha, 
  const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trmm");
#endif
    if( side == Left && shape == Lower )
    {
        if( orientation == Normal )
            BLAS::Internal::TrmmLLN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrmmLLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Left && shape == Upper )
    {
        if( orientation == Normal )
            BLAS::Internal::TrmmLUN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrmmLUT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Lower )
    {
        if( orientation == Normal )
            BLAS::Internal::TrmmRLN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrmmRLT( orientation, diagonal, alpha, A, X );
    }
    else if( side == Right && shape == Upper )
    {
        if( orientation == Normal )
            BLAS::Internal::TrmmRUN( diagonal, alpha, A, X );
        else
            BLAS::Internal::TrmmRUT( orientation, diagonal, alpha, A, X );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& A,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& A,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& A,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape, 
  const Orientation orientation, 
  const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& A,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

