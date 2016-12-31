/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

namespace El {

template<typename T>
void Ger
( T alpha,
  const Matrix<T>& x,
  const Matrix<T>& y,
        Matrix<T>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( ( x.Height() != 1 && x.Width() != 1 ) ||
          ( y.Height() != 1 && y.Width() != 1 ) )
          LogicError("x and y must be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
      if( xLength != A.Height() || yLength != A.Width() )
          LogicError
          ("Nonconformal Ger:\n",
           DimsString(x,"x"),"\n",DimsString(y,"y"),"\n",DimsString(A,"A"));
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
                   A.Buffer(), A.LDim() );
}

template<typename T>
void Ger
( T alpha,
  const AbstractDistMatrix<T>& x,
  const AbstractDistMatrix<T>& y,
        AbstractDistMatrix<T>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( APre, x, y );
      if( ( x.Width() != 1 && x.Height() != 1 ) ||
          ( y.Width() != 1 && y.Height() != 1 )   )
          LogicError("x and y are assumed to be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
      if( APre.Height() != xLength || APre.Width() != yLength )
          LogicError
          ("Nonconformal Ger:\n",
           DimsString(APre,"A"),"\n",DimsString(x,"x"),"\n",
           DimsString(y,"y"));
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> y_MR_STAR(g);
        y_MR_STAR.AlignWith( A );
        y_MR_STAR = y;

        Ger
        ( alpha, x_MC_STAR.LockedMatrix(),
                 y_MR_STAR.LockedMatrix(),
                 A.Matrix() );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,STAR,MR> y_STAR_MR(g);
        y_STAR_MR.AlignWith( A );
        y_STAR_MR = y;

        Ger
        ( alpha, x_MC_STAR.LockedMatrix(),
                 y_STAR_MR.LockedMatrix(),
                 A.Matrix() );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x;

        DistMatrix<T,MR,STAR> y_MR_STAR(g);
        y_MR_STAR.AlignWith( A );
        y_MR_STAR = y;

        Ger
        ( alpha, x_STAR_MC.LockedMatrix(),
                 y_MR_STAR.LockedMatrix(),
                 A.Matrix() );
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x;

        DistMatrix<T,STAR,MR> y_STAR_MR(g);
        y_STAR_MR.AlignWith( A );
        y_STAR_MR = y;

        Ger
        ( alpha, x_STAR_MC.LockedMatrix(),
                 y_STAR_MR.LockedMatrix(),
                 A.Matrix() );
    }
}

template<typename T>
void LocalGer
( T alpha,
  const AbstractDistMatrix<T>& x,
  const AbstractDistMatrix<T>& y,
        AbstractDistMatrix<T>& A )
{
    EL_DEBUG_CSE
    // TODO(poulson): Add error checking here
    Ger( alpha, x.LockedMatrix(), y.LockedMatrix(), A.Matrix() );
}

#define PROTO(T) \
  template void Ger \
  ( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A ); \
  template void Ger \
  ( T alpha, const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y, \
                   AbstractDistMatrix<T>& A ); \
  template void LocalGer \
  ( T alpha, const AbstractDistMatrix<T>& x, const AbstractDistMatrix<T>& y, \
                   AbstractDistMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
