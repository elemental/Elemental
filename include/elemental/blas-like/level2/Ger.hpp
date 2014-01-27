/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_GER_HPP
#define ELEM_GER_HPP

namespace elem {

template<typename T,Dist xColDist,Dist xRowDist,
                    Dist yColDist,Dist yRowDist,
                    Dist AColDist,Dist ARowDist>
inline void LocalGer
( T alpha, const DistMatrix<T,xColDist,xRowDist>& x,
           const DistMatrix<T,yColDist,yRowDist>& y,
                 DistMatrix<T,AColDist,ARowDist>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGer"))
    // TODO: Add error checking here
    Ger( alpha, x.LockedMatrix(), y.LockedMatrix(), A.Matrix() );
}

template<typename T>
inline void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("Ger");
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
inline void
Gerc( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{ Ger( alpha, x, y, A ); }

template<typename T>
inline void
Ger
( T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{
    DEBUG_ONLY(
        CallStackEntry cse("Ger");
        if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
            LogicError("{A,x,y} must be distributed over the same grid");
        if( ( x.Width() != 1 && x.Height() != 1 ) ||
            ( y.Width() != 1 && y.Height() != 1 )   )
            LogicError("x and y are assumed to be vectors");
        const Int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
        const Int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
        if( A.Height() != xLength || A.Width() != yLength )
            LogicError
            ("Nonconformal Ger:\n",
             DimsString(A,"A"),"\n",DimsString(x,"x"),"\n",DimsString(y,"y"));
    )
    const Grid& g = A.Grid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> y_MR_STAR(g);

        // Begin the algoritm
        x_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        y_MR_STAR = y;
        Ger
        ( alpha, x_MC_STAR.LockedMatrix(),
                 y_MR_STAR.LockedMatrix(),
                 A.Matrix() );
        //--------------------------------------------------------------------//
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,  STAR> x_MC_STAR(g);
        DistMatrix<T,STAR,MR  > y_STAR_MR(g);

        // Begin the algorithm
        x_MC_STAR.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        y_STAR_MR = y;
        Ger
        ( alpha, x_MC_STAR.LockedMatrix(),
                 y_STAR_MR.LockedMatrix(),
                 A.Matrix() );
        //--------------------------------------------------------------------//
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC  > x_STAR_MC(g);
        DistMatrix<T,MR,  STAR> y_MR_STAR(g);

        // Begin the algorithm
        x_STAR_MC.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        y_MR_STAR = y;
        Ger
        ( alpha, x_STAR_MC.LockedMatrix(),
                 y_MR_STAR.LockedMatrix(),
                 A.Matrix() );
        //--------------------------------------------------------------------//
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> y_STAR_MR(g);

        // Begin the algorithm
        x_STAR_MC.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MC = x;
        y_STAR_MR = y;
        Ger
        ( alpha, x_STAR_MC.LockedMatrix(),
                 y_STAR_MR.LockedMatrix(),
                 A.Matrix() );
        //--------------------------------------------------------------------//
    }
}

template<typename T>
inline void
Gerc
( T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{ Ger( alpha, x, y, A ); }

} // namespace elem

#endif // ifndef ELEM_GER_HPP
