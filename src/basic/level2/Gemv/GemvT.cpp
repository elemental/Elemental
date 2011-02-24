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
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

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
elemental::basic::internal::GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::internal::GemvT");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw logic_error( "{A,x,y} must be distributed over the same grid." );
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw logic_error( "GemvT expects x and y to be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Width() != yLength )
    {
        ostringstream msg;
        msg << "Nonconformal GemvT: " << endl
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
        DistMatrix<T,MR,Star> z_MR_Star(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T,MC,MR  > z(g);

        // Start the algorithm
        basic::Scal( beta, y );
        x_MC_Star.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        basic::Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(), 
                 x_MC_Star.LockedLocalMatrix(),
          (T)0,  z_MR_Star.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_Star );
        z = z_MR_MC;
        basic::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z_MR_MC.FreeAlignments();
        z.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(g);
        DistMatrix<T,MR,Star> z_MR_Star(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T,MC,MR  > zTrans(g);

        // Start the algorithm
        basic::Scal( beta, y );
        x_MC_Star.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        basic::Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(),
                 x_MC_Star.LockedLocalMatrix(),
          (T)0,  z_MR_Star.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_Star );
        basic::Trans( z_MR_MC, zTrans );
        basic::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z_MR_MC.FreeAlignments();
        zTrans.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(g);
        DistMatrix<T,MR,  Star> z_MR_Star(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);
        DistMatrix<T,MC,  MR  > z(g);

        // Start the algorithm
        basic::Scal( beta, y );
        x_Star_MC.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        z.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        basic::Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(), 
                 x_Star_MC.LockedLocalMatrix(),
          (T)0,  z_MR_Star.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_Star );
        z = z_MR_MC;
        basic::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z_MR_MC.FreeAlignments();
        z.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC  > x_Star_MC(g);
        DistMatrix<T,MR,  Star> z_MR_Star(g);
        DistMatrix<T,MR,  MC  > z_MR_MC(g);
        DistMatrix<T,MC,  MR  > zTrans(g);

        // Start the algorithm
        basic::Scal( beta, y );
        x_Star_MC.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z_MR_Star.ResizeTo( A.Width(), 1 );
        z_MR_MC.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_Star_MC = x;
        basic::Gemv
        ( orientation,
          alpha, A.LockedLocalMatrix(),
                 x_Star_MC.LockedLocalMatrix(),
          (T)0,  z_MR_Star.LocalMatrix() );
        z_MR_MC.SumScatterFrom( z_MR_Star );
        basic::Trans( z_MR_MC, zTrans );
        basic::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z_MR_MC.FreeAlignments();
        zTrans.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::GemvT
( Orientation orientation,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& x,
  float beta,        DistMatrix<float,MC,MR>& y );

template void elemental::basic::internal::GemvT
( Orientation orientation,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& x,
  double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::GemvT
( Orientation orientation,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& x,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void elemental::basic::internal::GemvT
( Orientation orientation,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& x,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif

