/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GEMV_T_HPP
#define ELEM_GEMV_T_HPP

#include ELEM_AXPY_INC
#include ELEM_SCALE_INC
#include ELEM_TRANSPOSE_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace internal {

template<typename T>
inline void
GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::GemvT");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( ( x.Width() != 1 && x.Height() != 1 ) ||
            ( y.Width() != 1 && y.Height() != 1 )   )
            LogicError("GemvT expects x and y to be vectors");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( A.Height() != xLength || A.Width() != yLength )
            LogicError
            ("Nonconformal GemvT: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
    )
    const Grid& g = A.Grid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);

        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_MC.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Axpy( T(1), z_MR_MC, y );
        //--------------------------------------------------------------------//
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> zTrans(g);

        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC  > x_STAR_MC(g);
        DistMatrix<T,MR,  STAR> z_MR_STAR(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);

        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_MC.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_STAR_MC, T(0), z_MR_STAR );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Axpy( T(1), z_MR_MC, y );
        //--------------------------------------------------------------------//
    }
    else
    {
        DistMatrix<T,STAR,MC  > x_STAR_MC(g);
        DistMatrix<T,MR,  STAR> z_MR_STAR(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);
        DistMatrix<T> zTrans(g);

        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        Zeros( z_MR_STAR, A.Width(), 1 );
        LocalGemv( orientation, alpha, A, x_STAR_MC, T(0), z_MR_STAR );
        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
}

template<typename T>
inline void
GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::GemvT");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( x.Width() != 1 || y.Width() != 1 )
            LogicError("GemvT expects x and y to be column vectors");
        if( A.Height() != x.Height() || A.Width() != y.Height() )
            LogicError
            ("Nonconformal GemvT: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
    )
    const Grid& g = A.Grid();

    DistMatrix<T,MC,STAR> x_MC_STAR(g);
    DistMatrix<T,MR,STAR> z_MR_STAR(g);
    DistMatrix<T,VR,STAR> z_VR_STAR(g);

    Scale( beta, y );
    x_MC_STAR.AlignWith( A );
    z_MR_STAR.AlignWith( A );
    z_VR_STAR.AlignWith( A );
    //--------------------------------------------------------------------//
    x_MC_STAR = x;
    Zeros( z_MR_STAR, A.Width(), 1 );
    LocalGemv( orientation, alpha, A, x_MC_STAR, T(0), z_MR_STAR );
    z_VR_STAR.PartialColSumScatterFrom( z_MR_STAR );
    Axpy( T(1), z_VR_STAR, y );
    //--------------------------------------------------------------------//
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_GEMV_T_HPP
