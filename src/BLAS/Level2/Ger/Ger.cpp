/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Ger
( const T alpha, const DistMatrix<T,MC,MR>& x,
                 const DistMatrix<T,MC,MR>& y,
                       DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::Ger");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw "{A,x,y} must be distributed over the same grid.";
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
    {
        throw "x and y are assumed to be vectors.";
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        ostringstream msg;
        msg << "Nonconformal Ger: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        const string s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = A.GetGrid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> y_MR_Star(grid);

        // Begin the algoritm
        x_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_MR_Star = y;
        BLAS::Ger( alpha, x_MC_Star.LockedLocalMatrix(),
                          y_MR_Star.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,  Star> x_MC_Star(grid);
        DistMatrix<T,Star,MR  > y_Star_MR(grid);

        // Begin the algorithm
        x_MC_Star.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_Star_MR = y;
        BLAS::Ger( alpha, x_MC_Star.LockedLocalMatrix(),
                          y_Star_MR.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        y_Star_MR.FreeConstraints();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(grid);
        DistMatrix<T,MR,  Star> y_MR_Star(grid);

        // Begin the algorithm
        x_Star_MC.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        y_MR_Star = y;
        BLAS::Ger( alpha, x_Star_MC.LockedLocalMatrix(),
                          y_MR_Star.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(grid);
        DistMatrix<T,Star,MR> y_Star_MR(grid);

        // Begin the algorithm
        x_Star_MC.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        y_Star_MR = y;
        BLAS::Ger( alpha, x_Star_MC.LockedLocalMatrix(),
                          y_Star_MR.LockedLocalMatrix(),
                          A.LocalMatrix()               );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        y_Star_MR.FreeConstraints();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Ger
( const float alpha, const DistMatrix<float,MC,MR>& x,
                     const DistMatrix<float,MC,MR>& y,
                           DistMatrix<float,MC,MR>& A );

template void Elemental::BLAS::Ger
( const double alpha, const DistMatrix<double,MC,MR>& x,
                      const DistMatrix<double,MC,MR>& y,
                            DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Ger
( const scomplex alpha, const DistMatrix<scomplex,MC,MR>& x,
                        const DistMatrix<scomplex,MC,MR>& y,
                              DistMatrix<scomplex,MC,MR>& A );

template void Elemental::BLAS::Ger
( const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& x,
                        const DistMatrix<dcomplex,MC,MR>& y,
                              DistMatrix<dcomplex,MC,MR>& A );
#endif
