/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

#include "./Symv/L.hpp"
#include "./Symv/U.hpp"

namespace El {

template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha,
  const Matrix<T>& A,
  const Matrix<T>& x,
  T beta,
        Matrix<T>& y,
  bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( ( x.Height() != 1 && x.Width() != 1 ) ||
          ( y.Height() != 1 && y.Width() != 1 ) )
          LogicError("x and y must be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
      if( A.Height() != xLength || A.Height() != yLength )
          LogicError("A must conform with x and y");
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const Int m = A.Height();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( conjugate )
    {
        blas::Hemv
        ( uploChar, m,
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
          beta,  y.Buffer(), incy );
    }
    else
    {
        blas::Symv
        ( uploChar, m,
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
          beta,  y.Buffer(), incy );
    }
}

template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& x,
  T beta,
        AbstractDistMatrix<T>& yPre,
  bool conjugate,
  const SymvCtrl<T>& ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, x, yPre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      if( ( x.Width() != 1 && x.Height() != 1 ) ||
          ( yPre.Width() != 1 && yPre.Height() != 1 ) )
          LogicError("x and y are assumed to be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( yPre.Width()==1 ? yPre.Height() : yPre.Width() );
      if( APre.Height() != xLength || APre.Height() != yLength )
          LogicError
          ("Nonconformal Symv: \n",DimsString(APre,"A"),"\n",
           DimsString(x,"x"),"\n",DimsString(yPre,"y"));
    )
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<T,T,MC,MR> yProx( yPre );
    auto& A = AProx.GetLocked();
    auto& y = yProx.Get();

    y *= beta;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x_MC_STAR;

        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MC_STAR.Resize( y.Height(), 1 );
        z_MR_STAR.Resize( y.Height(), 1 );
        Zero( z_MC_STAR );
        Zero( z_MR_STAR );
        if( uplo == LOWER )
        {
            symv::LocalColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
              ctrl );
        }
        else
        {
            symv::LocalColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
              ctrl );
        }

        DistMatrix<T,MR,MC> z_MR_MC(g);
        Contract( z_MR_STAR, z_MR_MC );

        DistMatrix<T> z(g);
        z.AlignWith( y );
        z = z_MR_MC;
        AxpyContract( T(1), z_MC_STAR, z );
        Axpy( T(1), z, y );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x_MC_STAR;

        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MC_STAR.Resize( y.Width(), 1 );
        z_MR_STAR.Resize( y.Width(), 1 );
        Zero( z_MC_STAR );
        Zero( z_MR_STAR );
        if( uplo == LOWER )
        {
            symv::LocalColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
              ctrl );
        }
        else
        {
            symv::LocalColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
              ctrl );
        }

        DistMatrix<T> z(g);
        z.AlignWith( y );
        Contract( z_MC_STAR, z );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC = z;
        AxpyContract( T(1), z_MR_STAR, z_MR_MC );

        DistMatrix<T> zTrans(g);
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;

        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x_STAR_MR;

        DistMatrix<T,STAR,MC> z_STAR_MC(g);
        DistMatrix<T,STAR,MR> z_STAR_MR(g);
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.AlignWith( A );
        z_STAR_MC.Resize( 1, y.Height() );
        z_STAR_MR.Resize( 1, y.Height() );
        Zero( z_STAR_MC );
        Zero( z_STAR_MR );
        if( uplo == LOWER )
        {
            symv::LocalRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
              ctrl );
        }
        else
        {
            symv::LocalRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
              ctrl );
        }

        DistMatrix<T> z(g);
        z.AlignWith( y );
        Contract( z_STAR_MR, z );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC = z;
        AxpyContract( T(1), z_STAR_MC, z_MR_MC );

        DistMatrix<T> zTrans(g);
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
    }
    else
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;

        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x_STAR_MR;

        DistMatrix<T,STAR,MC> z_STAR_MC(g);
        DistMatrix<T,STAR,MR> z_STAR_MR(g);
        z_STAR_MR.AlignWith( A );
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.Resize( 1, y.Width() );
        z_STAR_MC.Resize( 1, y.Width() );
        Zero( z_STAR_MR );
        Zero( z_STAR_MC );
        if( uplo == LOWER )
        {
            symv::LocalRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
              ctrl );
        }
        else
        {
            symv::LocalRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
              ctrl );
        }

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        Contract( z_STAR_MC, z_MR_MC );

        DistMatrix<T> z(g);
        z.AlignWith( y );
        z = z_MR_MC;
        AxpyContract( T(1), z_STAR_MR, z );
        Axpy( T(1), z, y );
    }
}

namespace symv {

template<typename T>
void LocalColAccumulate
( UpperOrLower uplo, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR, bool conjugate,
  const SymvCtrl<T>& ctrl )
{
    if( uplo == LOWER )
        LocalColAccumulateL
        ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
          ctrl );
    else
        LocalColAccumulateU
        ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
          ctrl );
}

template<typename T>
void LocalRowAccumulate
( UpperOrLower uplo, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR, bool conjugate,
  const SymvCtrl<T>& ctrl )
{
    if( uplo == LOWER )
        LocalRowAccumulateL
        ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
          ctrl );
    else
        LocalRowAccumulateU
        ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate,
          ctrl );
}

} // namespace symv

#define PROTO(T) \
  template void Symv \
  ( UpperOrLower uplo, \
    T alpha, \
    const Matrix<T>& A, \
    const Matrix<T>& x, \
    T beta, \
         Matrix<T>& y, \
    bool conjugate ); \
  template void Symv \
  ( UpperOrLower uplo, \
    T alpha, \
    const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& x, \
    T beta, \
          AbstractDistMatrix<T>& y, \
    bool conjugate, \
    const SymvCtrl<T>& ctrl ); \
  template void symv::LocalColAccumulate \
  ( UpperOrLower uplo, \
    T alpha, \
    const DistMatrix<T>& A, \
    const DistMatrix<T,MC,STAR>& x_MC_STAR, \
    const DistMatrix<T,MR,STAR>& x_MR_STAR, \
          DistMatrix<T,MC,STAR>& z_MC_STAR, \
          DistMatrix<T,MR,STAR>& z_MR_STAR, bool conjugate, \
    const SymvCtrl<T>& ctrl ); \
  template void symv::LocalRowAccumulate \
  ( UpperOrLower uplo, \
    T alpha, \
    const DistMatrix<T>& A, \
    const DistMatrix<T,STAR,MC>& x_STAR_MC, \
    const DistMatrix<T,STAR,MR>& x_STAR_MR, \
          DistMatrix<T,STAR,MC>& z_STAR_MC, \
          DistMatrix<T,STAR,MR>& z_STAR_MR, bool conjugate, \
    const SymvCtrl<T>& ctrl );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
