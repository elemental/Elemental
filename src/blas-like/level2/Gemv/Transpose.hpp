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
inline void
Transpose
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemv::Transpose");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( ( x.Width() != 1 && x.Height() != 1 ) ||
            ( y.Width() != 1 && y.Height() != 1 )   )
            LogicError("Expected x and y to be vectors");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( A.Height() != xLength || A.Width() != yLength )
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
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MR_STAR.AlignWith( A );
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Axpy( T(1), z_MR_MC, y );
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        x_MC_STAR.AlignWith( A );
        x_MC_STAR = x;

        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MR_STAR.AlignWith( A );
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );

        DistMatrix<T> zTrans(g);
        zTrans.AlignWith( y );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x;

        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MR_STAR.AlignWith( A );
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_STAR_MC, T(0), z_MR_STAR );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Axpy( T(1), z_MR_MC, y );
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        x_STAR_MC.AlignWith( A );
        x_STAR_MC = x;

        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        z_MR_STAR.AlignWith( A );
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_STAR_MC, T(0), z_MR_STAR );

        DistMatrix<T,MR,MC> z_MR_MC(g);
        z_MR_MC.AlignWith( y );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );

        DistMatrix<T> zTrans(g);
        zTrans.AlignWith( y );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
    }
}

template<typename T>
inline void
Transpose
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("gemv::Transpose");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("Expected x and y to be column vectors");
        if( A.Height() != x.Height() || A.Width() != y.Height() )
            LogicError
            ("Nonconformal: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
    )
    const Grid& g = A.Grid();
    Scale( beta, y );

    DistMatrix<T,MC,STAR> x_MC_STAR(g);
    x_MC_STAR.AlignWith( A );
    x_MC_STAR = x;

    DistMatrix<T,MR,STAR> z_MR_STAR(g);
    z_MR_STAR.AlignWith( A );
    Zeros( z_MR_STAR, A.Width(), 1 );
    LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );

    DistMatrix<T,VR,STAR> z_VR_STAR(g);
    z_VR_STAR.AlignWith( A );
    z_VR_STAR.PartialColSumScatterFrom( z_MR_STAR );
    Axpy( T(1), z_VR_STAR, y );
}

} // namespace gemv
} // namespace El
