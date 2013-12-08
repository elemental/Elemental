/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_GEMV_N_HPP
#define ELEM_BLAS_GEMV_N_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"

#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace internal {

template<typename T>
inline void
GemvN
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    CallStackEntry cse("internal::GemvN");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        LogicError("{A,x,y} must be distributed over the same grid");
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        LogicError("x and y are assumed to be vectors");
    const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != yLength || A.Width() != xLength )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemvN: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        LogicError( msg.str() );
    }
#endif
    const Grid& g = A.Grid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);

        Scale( beta, y );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MR_STAR = x;
        Zeros( z_MC_STAR, A.Height(), 1 );
        LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
        y.SumScatterUpdate( T(1), z_MC_STAR );
        //--------------------------------------------------------------------//
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T> z(g), zTrans(g);

        Scale( beta, y );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MR_STAR = x;
        Zeros( z_MC_STAR, A.Height(), 1 );
        LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
        z.SumScatterFrom( z_MC_STAR );
        Transpose( z, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MR  > x_STAR_MR(g);
        DistMatrix<T,MC,  STAR> z_MC_STAR(g);

        Scale( beta, y );
        x_STAR_MR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        Zeros( z_MC_STAR, A.Height(), 1 );
        LocalGemv( NORMAL, alpha, A, x_STAR_MR, T(0), z_MC_STAR );
        y.SumScatterUpdate( T(1), z_MC_STAR );
        //--------------------------------------------------------------------//
    }
    else
    {
        DistMatrix<T,STAR,MR  > x_STAR_MR(g);
        DistMatrix<T,MC,  STAR> z_MC_STAR(g);
        DistMatrix<T> z(g), zTrans(g);

        Scale( beta, y );
        x_STAR_MR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        Zeros( z_MC_STAR, A.Height(), 1 );
        LocalGemv( NORMAL, alpha, A, x_STAR_MR, T(0), z_MC_STAR );
        z.SumScatterFrom( z_MC_STAR );
        Transpose( z, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
}

template<typename T>
inline void
GemvN
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y )
{
#ifndef RELEASE
    CallStackEntry cse("internal::GemvN");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        LogicError("{A,x,y} must be distributed over the same grid");
    if( x.Width() != 1 || y.Width() != 1 )
        LogicError("x and y are assumed to be column vectors");
    if( A.Height() != y.Height() || A.Width() != x.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemvN: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        LogicError( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    DistMatrix<T,MR,STAR> x_MR_STAR(g);
    DistMatrix<T,MC,STAR> z_MC_STAR(g);

    Scale( beta, y );
    x_MR_STAR.AlignWith( A );
    z_MC_STAR.AlignWith( A );
    //--------------------------------------------------------------------//
    x_MR_STAR = x;
    Zeros( z_MC_STAR, A.Height(), 1 );
    LocalGemv( NORMAL, alpha, A, x_MR_STAR, T(0), z_MC_STAR );
    y.SumScatterUpdate( T(1), z_MC_STAR );
    //--------------------------------------------------------------------//
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_BLAS_GEMV_N_HPP
