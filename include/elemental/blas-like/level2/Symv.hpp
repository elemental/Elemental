/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SYMV_HPP
#define ELEM_SYMV_HPP

#include ELEM_AXPY_INC
#include ELEM_SCALE_INC
#include ELEM_TRANSPOSE_INC

#include ELEM_ZEROS_INC

#include "./Symv/L.hpp"
#include "./Symv/U.hpp"

namespace elem {

template<typename T>
inline void
Symv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Symv");
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
inline void
Symv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Symv");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( ( x.Width() != 1 && x.Height() != 1 ) ||
            ( y.Width() != 1 && y.Height() != 1 ) )
            LogicError("x and y are assumed to be vectors");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( A.Height() != xLength || A.Height() != yLength )
            LogicError
            ("Nonconformal Symv: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x ~ ",x.Height()," x ",x.Width(),"\n",
             "  y ~ ",y.Height()," x ",y.Width(),"\n");
    )
    const Grid& g = A.Grid();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,STAR> x_MC_STAR(g), z_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g), z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> z(g);

        // Begin the algoritm
        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z.AlignWith( y );
        Zeros( z_MC_STAR, y.Height(), 1 );
        Zeros( z_MR_STAR, y.Height(), 1 );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            internal::LocalSymvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate );
        }
        else
        {
            internal::LocalSymvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate );
        }

        z_MR_MC.RowSumScatterFrom( z_MR_STAR );
        z = z_MR_MC;
        z.RowSumScatterUpdate( T(1), z_MC_STAR );
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,STAR> x_MC_STAR(g), z_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g), z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> z(g), zTrans(g);

        // Begin the algoritm
        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        Zeros( z_MC_STAR, y.Width(), 1 );
        Zeros( z_MR_STAR, y.Width(), 1 );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            internal::LocalSymvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate );
        }
        else
        {
            internal::LocalSymvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate );
        }

        z.RowSumScatterFrom( z_MC_STAR );
        z_MR_MC = z;
        z_MR_MC.RowSumScatterUpdate( T(1), z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC> x_STAR_MC(g), z_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g), z_STAR_MR(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);
        DistMatrix<T> z(g), zTrans(g);

        // Begin the algoritm
        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        Zeros( z_STAR_MC, 1, y.Height() );
        Zeros( z_STAR_MR, 1, y.Height() );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            internal::LocalSymvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate );
        }
        else
        {
            internal::LocalSymvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate );
        }

        z.ColSumScatterFrom( z_STAR_MR );
        z_MR_MC = z;
        z_MR_MC.ColSumScatterUpdate( T(1), z_STAR_MC );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC> x_STAR_MC(g), z_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g), z_STAR_MR(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);
        DistMatrix<T> z(g);

        // Begin the algoritm
        Scale( beta, y );
        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        Zeros( z_STAR_MC, 1, y.Width() );
        Zeros( z_STAR_MR, 1, y.Width() );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            internal::LocalSymvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate );
        }
        else
        {
            internal::LocalSymvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR, conjugate );
        }

        z_MR_MC.ColSumScatterFrom( z_STAR_MC );
        z = z_MR_MC;
        z.ColSumScatterUpdate( T(1), z_STAR_MR );
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
    }
}

} // namespace elem

#endif // ifndef ELEM_SYMV_HPP
