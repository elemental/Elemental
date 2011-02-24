/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/basic.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::utilities;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

template<typename T>
void
elemental::basic::Her2
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("basic::Her2");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
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
    const Grid& g = A.Grid();

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int r = g.Height();
    const int c = g.Width();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();

    if( x.Width() == 1 && y.Width() == 1 )
    {
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(y_MR_Star.GetLocalEntry(jLoc,0)) + 
                                        y_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(x_MR_Star.GetLocalEntry(jLoc,0)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(y_MR_Star.GetLocalEntry(jLoc,0)) + 
                                        y_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(x_MR_Star.GetLocalEntry(jLoc,0)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(y_Star_MR.GetLocalEntry(0,jLoc)) + 
                                        y_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(x_MR_Star.GetLocalEntry(jLoc,0)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(y_Star_MR.GetLocalEntry(0,jLoc)) +
                                        y_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(x_MR_Star.GetLocalEntry(jLoc,0)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_Star_MC.GetLocalEntry(0,iLoc)* 
                        elemental::Conj(y_MR_Star.GetLocalEntry(jLoc,0)) + 
                                        y_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(x_Star_MR.GetLocalEntry(0,jLoc)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(y_MR_Star.GetLocalEntry(jLoc,0)) +
                                        y_MC_Star.GetLocalEntry(iLoc,0)*
                        elemental::Conj(x_Star_MR.GetLocalEntry(0,jLoc)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(y_Star_MR.GetLocalEntry(0,jLoc)) + 
                                        y_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(x_Star_MR.GetLocalEntry(0,jLoc)) ) );
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
                    const T value = A.GetLocalEntry(iLoc,jLoc);
                    A.SetLocalEntry
                    ( iLoc, jLoc,
                      value + alpha*
                      (                 x_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(y_Star_MR.GetLocalEntry(0,jLoc)) + 
                                        y_Star_MC.GetLocalEntry(0,iLoc)*
                        elemental::Conj(x_Star_MR.GetLocalEntry(0,jLoc)) ) );
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

template void elemental::basic::Her2
( Shape shape, 
  float alpha, const DistMatrix<float,MC,MR>& x,
               const DistMatrix<float,MC,MR>& y,
                     DistMatrix<float,MC,MR>& A );

template void elemental::basic::Her2
( Shape shape,
  double alpha, const DistMatrix<double,MC,MR>& x,
                const DistMatrix<double,MC,MR>& y,
                      DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::Her2
( Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& x,
                  const DistMatrix<scomplex,MC,MR>& y,
                        DistMatrix<scomplex,MC,MR>& A );

template void elemental::basic::Her2
( Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& x,
                  const DistMatrix<dcomplex,MC,MR>& y,
                        DistMatrix<dcomplex,MC,MR>& A );
#endif

