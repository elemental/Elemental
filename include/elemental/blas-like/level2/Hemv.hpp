/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Hemv/L.hpp"
#include "./Hemv/U.hpp"

namespace elem {

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
        throw std::logic_error("A must conform with x and y");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 ) )
        throw std::logic_error("x and y are assumed to be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        std::ostringstream msg;
        msg << "Nonconformal Hemv: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T> z(g);

        // Begin the algoritm
        Scale( beta, y );
        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z.AlignWith( y );
        z_MC_STAR.ResizeTo( y.Height(), 1 );
        z_MR_STAR.ResizeTo( y.Height(), 1 );
        Zero( z_MC_STAR );
        Zero( z_MR_STAR );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }
        else
        {
            internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }

        z_MR_MC.SumScatterFrom( z_MR_STAR );
        z = z_MR_MC;
        z.SumScatterUpdate( T(1), z_MC_STAR );
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T,MR,STAR> z_MR_STAR(g);
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
        z_MC_STAR.ResizeTo( y.Width(), 1 );
        z_MR_STAR.ResizeTo( y.Width(), 1 );
        Zero( z_MC_STAR );
        Zero( z_MR_STAR );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }
        else
        {
            internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }

        z.SumScatterFrom( z_MC_STAR );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( T(1), z_MR_STAR );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
        z_MR_STAR.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        DistMatrix<T,STAR,MC> z_STAR_MC(g);
        DistMatrix<T,STAR,MR> z_STAR_MR(g);
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
        z_STAR_MC.ResizeTo( 1, y.Height() );
        z_STAR_MR.ResizeTo( 1, y.Height() );
        Zero( z_STAR_MC );
        Zero( z_STAR_MR );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            internal::LocalHemvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }
        else
        {
            internal::LocalHemvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }

        z.SumScatterFrom( z_STAR_MR );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( T(1), z_STAR_MC );
        Transpose( z_MR_MC, zTrans );
        Axpy( T(1), zTrans, y );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
        z_STAR_MC.FreeAlignments();
        z_STAR_MR.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        DistMatrix<T,STAR,MC> z_STAR_MC(g);
        DistMatrix<T,STAR,MR> z_STAR_MR(g);
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
        z_STAR_MC.ResizeTo( 1, y.Width() );
        z_STAR_MR.ResizeTo( 1, y.Width() );
        Zero( z_STAR_MC );
        Zero( z_STAR_MR );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            internal::LocalHemvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }
        else
        {
            internal::LocalHemvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }

        z_MR_MC.SumScatterFrom( z_STAR_MC );
        z = z_MR_MC;
        z.SumScatterUpdate( T(1), z_STAR_MR );
        Axpy( T(1), z, y );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
        z_STAR_MC.FreeAlignments();
        z_STAR_MR.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
