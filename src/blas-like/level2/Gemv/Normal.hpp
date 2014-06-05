/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace gemv {

template<typename T>
inline void Normal
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemv::Normal");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( ( x.Width() != 1 && x.Height() != 1 ) ||
            ( y.Width() != 1 && y.Height() != 1 )   )
            LogicError("x and y are assumed to be vectors");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( A.Height() != yLength || A.Width() != xLength )
            LogicError
            ("Nonconformal: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
    )
    const Grid& g = A.Grid();
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
        y.RowSumScatterUpdate( T(1), z_MC_STAR );
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
        z.RowSumScatterFrom( z_MC_STAR );
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
        y.RowSumScatterUpdate( T(1), z_MC_STAR );
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
        z.RowSumScatterFrom( z_MC_STAR );
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
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("x and y are assumed to be column vectors");
        if( A.Height() != y.Height() || A.Width() != x.Height() )
            LogicError
            ("Nonconformal: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
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
    y.PartialColSumScatterUpdate( T(1), z_MC_STAR );
}

} // namespace gemv
} // namespace El
