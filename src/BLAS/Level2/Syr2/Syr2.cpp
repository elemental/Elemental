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
#include "Elemental/BLAS.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::utilities;

template<typename T>
void
BLAS::Syr2
( const Shape shape,
  const T alpha, const DistMatrix<T,MC,MR>& x,
                 const DistMatrix<T,MC,MR>& y,
                       DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syr2");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw "{A,x,y} must be distributed over the same grid.";
    if( A.Height() != A.Width() )
        throw "A must be square.";
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        ostringstream msg;
        msg << "A must conform with x: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        throw msg.str();
    }
#endif
    const Grid& grid = A.GetGrid();

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int r = grid.Height();
    const int c = grid.Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,MC,Star> y_MC_Star(grid);
        DistMatrix<T,MR,Star> y_MR_Star(grid);

        x_MC_Star.AlignWith( A );
        x_MR_Star.AlignWith( A );
        y_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        x_MR_Star = x_MC_Star;
        y_MC_Star = y;
        y_MR_Star = y_MC_Star;

        if( shape == Lower )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);
                for( int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_MC_Star.LocalEntry(iLoc,0) *
                                  y_MR_Star.LocalEntry(jLoc,0) +
                                  y_MC_Star.LocalEntry(iLoc,0) *
                                  x_MR_Star.LocalEntry(jLoc,0)   );
                }
            }
        }
        else
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);
                for( int iLoc=0; iLoc<heightToDiag; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_MC_Star.LocalEntry(iLoc,0) *
                                  y_MR_Star.LocalEntry(jLoc,0) +
                                  y_MC_Star.LocalEntry(iLoc,0) *
                                  x_MR_Star.LocalEntry(jLoc,0)   );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        x_MR_Star.FreeConstraints();
        y_MC_Star.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(grid);
        DistMatrix<T,MR,Star> x_MR_Star(grid);
        DistMatrix<T,Star,MC> y_Star_MC(grid);
        DistMatrix<T,Star,MR> y_Star_MR(grid);

        x_MC_Star.AlignWith( A );
        x_MR_Star.AlignWith( A );
        y_Star_MC.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        x_MR_Star = x_MC_Star;
        y_Star_MR = y;
        y_Star_MC = y_Star_MR;

        if( shape == Lower )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);
                for( int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_MC_Star.LocalEntry(iLoc,0) *
                                  y_Star_MR.LocalEntry(0,jLoc) +
                                  y_Star_MC.LocalEntry(0,iLoc) *
                                  x_MR_Star.LocalEntry(jLoc,0)   );
                }
            }
        }
        else
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);
                for( int iLoc=0; iLoc<heightToDiag; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_MC_Star.LocalEntry(iLoc,0) *
                                  y_Star_MR.LocalEntry(0,jLoc) +
                                  y_Star_MC.LocalEntry(0,iLoc) *
                                  x_MR_Star.LocalEntry(jLoc,0)   );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_MC_Star.FreeConstraints();
        x_MR_Star.FreeConstraints();
        y_Star_MC.FreeConstraints();
        y_Star_MR.FreeConstraints();

    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(grid);
        DistMatrix<T,Star,MR> x_Star_MR(grid);
        DistMatrix<T,MC,Star> y_MC_Star(grid);
        DistMatrix<T,MR,Star> y_MR_Star(grid);

        x_Star_MC.AlignWith( A );
        x_Star_MR.AlignWith( A );
        y_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        x_Star_MC = x_Star_MR;
        y_MC_Star = y;
        y_MR_Star = y_MC_Star;

        if( shape == Lower )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);
                for( int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_Star_MC.LocalEntry(0,iLoc) *
                                  y_MR_Star.LocalEntry(jLoc,0) + 
                                  y_MC_Star.LocalEntry(iLoc,0) *
                                  x_Star_MR.LocalEntry(0,jLoc)   );
                }
            }
        }
        else
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);
                for( int iLoc=0; iLoc<heightToDiag; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_Star_MC.LocalEntry(0,iLoc) *
                                  y_MR_Star.LocalEntry(jLoc,0) +
                                  y_MC_Star.LocalEntry(iLoc,0) *
                                  x_Star_MR.LocalEntry(0,jLoc)   );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        x_Star_MR.FreeConstraints();
        y_MC_Star.FreeConstraints();
        y_MR_Star.FreeConstraints();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(grid);
        DistMatrix<T,Star,MR> x_Star_MR(grid);
        DistMatrix<T,Star,MC> y_Star_MC(grid);
        DistMatrix<T,Star,MR> y_Star_MR(grid);

        x_Star_MC.AlignWith( A );
        x_Star_MR.AlignWith( A );
        y_Star_MC.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        x_Star_MC = x_Star_MR;
        y_Star_MR = y;
        y_Star_MC = y_Star_MR;

        if( shape == Lower )
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);
                for( int iLoc=heightAboveDiag; iLoc<localHeight; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_Star_MC.LocalEntry(0,iLoc) *
                                  y_Star_MR.LocalEntry(0,jLoc) + 
                                  y_Star_MC.LocalEntry(0,iLoc) *
                                  x_Star_MR.LocalEntry(0,jLoc)   );
                }
            }
        }
        else
        {
            for( int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const int j = rowShift + jLoc*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);
                for( int iLoc=0; iLoc<heightToDiag; ++iLoc )
                {
                    A.LocalEntry(iLoc,jLoc) += 
                        alpha * ( x_Star_MC.LocalEntry(0,iLoc) *
                                  y_Star_MR.LocalEntry(0,jLoc) +
                                  y_Star_MC.LocalEntry(0,iLoc) *
                                  x_Star_MR.LocalEntry(0,jLoc)   );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_Star_MC.FreeConstraints();
        x_Star_MR.FreeConstraints();
        y_Star_MC.FreeConstraints();
        y_Star_MR.FreeConstraints();
 
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Syr2
( const Shape shape, 
  const float alpha, const DistMatrix<float,MC,MR>& x,
                     const DistMatrix<float,MC,MR>& y,
                           DistMatrix<float,MC,MR>& A );

template void Elemental::BLAS::Syr2
( const Shape shape,
  const double alpha, const DistMatrix<double,MC,MR>& x,
                      const DistMatrix<double,MC,MR>& y,
                            DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Syr2
( const Shape shape,
  const scomplex alpha, const DistMatrix<scomplex,MC,MR>& x,
                        const DistMatrix<scomplex,MC,MR>& y,
                              DistMatrix<scomplex,MC,MR>& A );

template void Elemental::BLAS::Syr2
( const Shape shape,
  const dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& x,
                        const DistMatrix<dcomplex,MC,MR>& y,
                              DistMatrix<dcomplex,MC,MR>& A );
#endif
