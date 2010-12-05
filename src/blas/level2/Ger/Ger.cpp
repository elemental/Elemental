/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::Ger
( T alpha, const DistMatrix<T,MC,MR>& x,
           const DistMatrix<T,MC,MR>& y,
                 DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("blas::Internal::Ger");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw logic_error( "{A,x,y} must be distributed over the same grid." );
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw logic_error( "x and y are assumed to be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        ostringstream msg;
        msg << "Nonconformal Ger: " << endl
            << "  A ~ " << A.Height() << " x " << A.Width() << endl
            << "  x ~ " << x.Height() << " x " << x.Width() << endl
            << "  y ~ " << y.Height() << " x " << y.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = A.Grid();
    if( x.Width() == 1 && y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(g);
        DistMatrix<T,MR,Star> y_MR_Star(g);

        // Begin the algoritm
        x_MC_Star.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_MR_Star = y;
        blas::Ger
        ( alpha, x_MC_Star.LockedLocalMatrix(),
                 y_MR_Star.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        y_MR_Star.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,  Star> x_MC_Star(g);
        DistMatrix<T,Star,MR  > y_Star_MR(g);

        // Begin the algorithm
        x_MC_Star.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        y_Star_MR = y;
        blas::Ger
        ( alpha, x_MC_Star.LockedLocalMatrix(),
                 y_Star_MR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        y_Star_MR.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(g);
        DistMatrix<T,MR,  Star> y_MR_Star(g);

        // Begin the algorithm
        x_Star_MC.AlignWith( A );
        y_MR_Star.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        y_MR_Star = y;
        blas::Ger
        ( alpha, x_Star_MC.LockedLocalMatrix(),
                 y_MR_Star.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        y_MR_Star.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(g);
        DistMatrix<T,Star,MR> y_Star_MR(g);

        // Begin the algorithm
        x_Star_MC.AlignWith( A );
        y_Star_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        y_Star_MR = y;
        blas::Ger
        ( alpha, x_Star_MC.LockedLocalMatrix(),
                 y_Star_MR.LockedLocalMatrix(),
                 A.LocalMatrix() );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        y_Star_MR.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Ger
( float alpha, const DistMatrix<float,MC,MR>& x,
               const DistMatrix<float,MC,MR>& y,
                     DistMatrix<float,MC,MR>& A );

template void elemental::blas::Ger
( double alpha, const DistMatrix<double,MC,MR>& x,
                const DistMatrix<double,MC,MR>& y,
                      DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Ger
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& x,
                  const DistMatrix<scomplex,MC,MR>& y,
                        DistMatrix<scomplex,MC,MR>& A );

template void elemental::blas::Ger
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& x,
                  const DistMatrix<dcomplex,MC,MR>& y,
                        DistMatrix<dcomplex,MC,MR>& A );
#endif

