/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;

template<typename T>
void
elemental::blas::Syr2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Syr2");
    if( A.GetGrid() != x.GetGrid() || x.GetGrid() != y.GetGrid() )
        throw logic_error( "{A,x,y} must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        ostringstream msg;
        msg << "A must conform with x: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.GetGrid();

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(g);
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,MC,Star> y_MC_Star(g);
        DistMatrix<T,MR,Star> y_MR_Star(g);

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
        x_MC_Star.FreeAlignments();
        x_MR_Star.FreeAlignments();
        y_MC_Star.FreeAlignments();
        y_MR_Star.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(g);
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,Star,MC> y_Star_MC(g);
        DistMatrix<T,Star,MR> y_Star_MR(g);

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
        x_MC_Star.FreeAlignments();
        x_MR_Star.FreeAlignments();
        y_Star_MC.FreeAlignments();
        y_Star_MR.FreeAlignments();

    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(g);
        DistMatrix<T,Star,MR> x_Star_MR(g);
        DistMatrix<T,MC,Star> y_MC_Star(g);
        DistMatrix<T,MR,Star> y_MR_Star(g);

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
        x_Star_MC.FreeAlignments();
        x_Star_MR.FreeAlignments();
        y_MC_Star.FreeAlignments();
        y_MR_Star.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(g);
        DistMatrix<T,Star,MR> x_Star_MR(g);
        DistMatrix<T,Star,MC> y_Star_MC(g);
        DistMatrix<T,Star,MR> y_Star_MR(g);

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
        x_Star_MC.FreeAlignments();
        x_Star_MR.FreeAlignments();
        y_Star_MC.FreeAlignments();
        y_Star_MR.FreeAlignments();
 
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Syr2
( Shape shape, 
  float alpha, const DistMatrix<float,MC,MR>& x,
               const DistMatrix<float,MC,MR>& y,
                     DistMatrix<float,MC,MR>& A );

template void elemental::blas::Syr2
( Shape shape,
  double alpha, const DistMatrix<double,MC,MR>& x,
                const DistMatrix<double,MC,MR>& y,
                      DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Syr2
( Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& x,
                  const DistMatrix<scomplex,MC,MR>& y,
                        DistMatrix<scomplex,MC,MR>& A );

template void elemental::blas::Syr2
( Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& x,
                  const DistMatrix<dcomplex,MC,MR>& y,
                        DistMatrix<dcomplex,MC,MR>& A );
#endif

