/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemv {

template<typename T>
void Normal
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& x,
  T beta,
        AbstractDistMatrix<T>& yPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( APre, x, yPre );
      if( ( x.Width() != 1 && x.Height() != 1 ) ||
          ( yPre.Width() != 1 && yPre.Height() != 1 )   )
          LogicError("x and y are assumed to be vectors");
      const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
      const Int yLength = ( yPre.Width()==1 ? yPre.Height() : yPre.Width() );
      if( APre.Height() != yLength || APre.Width() != xLength )
          LogicError
          ("Nonconformal: \n",DimsString(APre,"A"),"\n",
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
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x;

        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.Resize( A.Height(), 1 );
        Zero( z_MC_STAR );
        LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
        AxpyContract( T(1), z_MC_STAR, y );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x;

        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.Resize( A.Height(), 1 );
        Zero( z_MC_STAR );
        LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );

        DistMatrix<T> z(g), zTrans(g);
        z.AlignWith( y );
        zTrans.AlignWith( y );
        Contract( z_MC_STAR, z );
        Transpose( z, zTrans );
        Axpy( T(1), zTrans, y );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;
        DistMatrix<T,MC,  STAR> z_MC_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.Resize( A.Height(), 1 );
        Zero( z_MC_STAR );
        LocalGemv( NORMAL, alpha, A, x_STAR_MR, T(0), z_MC_STAR );
        AxpyContract( T(1), z_MC_STAR, y );
    }
    else
    {
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        x_STAR_MR.AlignWith( A );
        x_STAR_MR = x;

        DistMatrix<T,MC,  STAR> z_MC_STAR(g);
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.Resize( A.Height(), 1 );
        Zero( z_MC_STAR );
        LocalGemv( NORMAL, alpha, A, x_STAR_MR, T(0), z_MC_STAR );

        DistMatrix<T> z(g), zTrans(g);
        z.AlignWith( y );
        zTrans.AlignWith( y );
        Contract( z_MC_STAR, z );
        Transpose( z, zTrans );
        Axpy( T(1), zTrans, y );
    }
}

template<typename T>
void Normal
( T alpha,
  const DistMatrix<T>& A,
  const AbstractDistMatrix<T>& x,
  T beta,
        DistMatrix<T,VC,STAR>& y )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, x, y );
      if( x.Width() != 1 || y.Width() != 1 )
          LogicError("x and y are assumed to be column vectors");
      if( A.Height() != y.Height() || A.Width() != x.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(A,"A"),"\n",
           DimsString(x,"x"),"\n",DimsString(y,"y"));
    )
    const Grid& g = A.Grid();
    y *= beta;

    DistMatrix<T,MR,STAR> x_MR_STAR(g);
    x_MR_STAR.AlignWith( A );
    x_MR_STAR = x;

    DistMatrix<T,MC,STAR> z_MC_STAR(g);
    z_MC_STAR.AlignWith( A );
    z_MC_STAR.Resize( A.Height(), 1 );
    Zero( z_MC_STAR );
    LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
    AxpyContract( T(1), z_MC_STAR, y );
}

} // namespace gemv
} // namespace El
