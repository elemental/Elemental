/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemv {

template<typename T>
inline void Normal
( T alpha, const AbstractDistMatrix<T>& APre,
           const AbstractDistMatrix<T>& xPre,
  T beta,        AbstractDistMatrix<T>& yPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemv::Normal");
        AssertSameGrids( APre, xPre, yPre );
        if( ( xPre.Width() != 1 && xPre.Height() != 1 ) ||
            ( yPre.Width() != 1 && yPre.Height() != 1 )   )
            LogicError("x and y are assumed to be vectors");
        const Int xLength = ( xPre.Width()==1 ? xPre.Height() : xPre.Width() );
        const Int yLength = ( yPre.Width()==1 ? yPre.Height() : yPre.Width() );
        if( APre.Height() != yLength || APre.Width() != xLength )
            LogicError
            ("Nonconformal: \n",DimsString(APre,"A"),"\n",
             DimsString(xPre,"x"),"\n",DimsString(yPre,"y"));
    )
    const Grid& g = APre.Grid();

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto xPtr = ReadProxy<T,MC,MR>( &xPre );      auto& x = *xPtr;
    auto yPtr = ReadWriteProxy<T,MC,MR>( &yPre ); auto& y = *yPtr;

    Scale( beta, y );
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        x_MR_STAR.AlignWith( A );
        x_MR_STAR = x;

        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        z_MC_STAR.AlignWith( A );
        Zeros( z_MC_STAR, A.Height(), 1 );
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
        Zeros( z_MC_STAR, A.Height(), 1 );
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
        Zeros( z_MC_STAR, A.Height(), 1 );
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
        Zeros( z_MC_STAR, A.Height(), 1 );
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
inline void Normal
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemv::Normal");
        AssertSameGrids( A, x, y );
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("x and y are assumed to be column vectors");
        if( A.Height() != y.Height() || A.Width() != x.Height() )
            LogicError
            ("Nonconformal: \n",DimsString(A,"A"),"\n",
             DimsString(x,"x"),"\n",DimsString(y,"y"));
    )
    const Grid& g = A.Grid();
    Scale( beta, y );

    DistMatrix<T,MR,STAR> x_MR_STAR(g);
    x_MR_STAR.AlignWith( A );
    x_MR_STAR = x;

    DistMatrix<T,MC,STAR> z_MC_STAR(g);
    z_MC_STAR.AlignWith( A );
    Zeros( z_MC_STAR, A.Height(), 1 );
    LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
    AxpyContract( T(1), z_MC_STAR, y );
}

} // namespace gemv
} // namespace El
