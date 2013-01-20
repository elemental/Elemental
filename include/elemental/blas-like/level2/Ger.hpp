/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_GER_HPP
#define BLAS_GER_HPP

namespace elem {

namespace internal {

template<typename T,Distribution xColDist,Distribution xRowDist,
                    Distribution yColDist,Distribution yRowDist,
                    Distribution AColDist,Distribution ARowDist>
inline void LocalGer
( T alpha, const DistMatrix<T,xColDist,xRowDist>& x,
           const DistMatrix<T,yColDist,yRowDist>& y,
                 DistMatrix<T,AColDist,ARowDist>& A )
{
#ifndef RELEASE
    PushCallStack("internal::LocalGer");
    // TODO: Add error checking here
#endif
    Ger( alpha, x.LockedLocalMatrix(), y.LockedLocalMatrix(), A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 ) )
        throw std::logic_error("x and y must be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger:\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n"
            << "  A ~ " << A.Height() << " x " << A.Width();
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    blas::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
                   A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
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
#ifndef RELEASE
    PushCallStack("Ger");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw std::logic_error("x and y are assumed to be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        std::ostringstream msg;
        msg << "Nonconformal Ger: \n"
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
        DistMatrix<T,MR,STAR> y_MR_STAR(g);

        // Begin the algoritm
        x_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        y_MR_STAR = y;
        Ger
        ( alpha, x_MC_STAR.LockedLocalMatrix(),
                 y_MR_STAR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        y_MR_STAR.FreeAlignments();
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
        ( alpha, x_MC_STAR.LockedLocalMatrix(),
                 y_STAR_MR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        y_STAR_MR.FreeAlignments();
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
        ( alpha, x_STAR_MC.LockedLocalMatrix(),
                 y_MR_STAR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        y_MR_STAR.FreeAlignments();
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
        ( alpha, x_STAR_MC.LockedLocalMatrix(),
                 y_STAR_MR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        y_STAR_MR.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Gerc
( T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{ Ger( alpha, x, y, A ); }

} // namespace elem

#endif // ifndef BLAS_GER_HPP
