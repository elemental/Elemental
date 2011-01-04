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
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

// Set up interface for managing tuning parameters
namespace {
int localHemvFloatBlocksize = 64;
int localHemvDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localHemvComplexFloatBlocksize = 64;
int localHemvComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX
} 

void elemental::blas::SetLocalHemvFloatBlocksize( int blocksize )
{ ::localHemvFloatBlocksize = blocksize; }

void elemental::blas::SetLocalHemvDoubleBlocksize( int blocksize )
{ ::localHemvDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
void elemental::blas::SetLocalHemvComplexFloatBlocksize( int blocksize )
{ ::localHemvComplexFloatBlocksize = blocksize; }

void elemental::blas::SetLocalHemvComplexDoubleBlocksize( int blocksize )
{ ::localHemvComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::blas::LocalHemvBlocksize<float>()
{ return ::localHemvFloatBlocksize; }

template<>
int
elemental::blas::LocalHemvBlocksize<double>()
{ return ::localHemvDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::blas::LocalHemvBlocksize<scomplex>()
{ return ::localHemvComplexFloatBlocksize; }

template<>
int
elemental::blas::LocalHemvBlocksize<dcomplex>()
{ return ::localHemvComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

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
elemental::blas::Hemv
( Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("blas::Hemv");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw logic_error( "{A,x,y} must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 ) )
        throw logic_error( "x and y are assumed to be vectors." );
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        ostringstream msg;
        msg << "Nonconformal Hemv: " << endl
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
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,MC,Star> z_MC_Star(g);
        DistMatrix<T,MR,Star> z_MR_Star(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T,MC,MR  > z(g);

        // Begin the algoritm
        blas::Scal( beta, y );
        x_MC_Star.AlignWith( A );
        x_MR_Star.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z.AlignWith( y );
        z_MC_Star.ResizeTo( y.Height(), 1 );
        z_MR_Star.ResizeTo( y.Height(), 1 );
        z_MC_Star.SetToZero();
        z_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        x_MR_Star = x_MC_Star;
        if( shape == Lower )
        {
            blas::internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_Star, x_MR_Star, z_MC_Star, z_MR_Star );
        }
        else
        {
            blas::internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_Star, x_MR_Star, z_MC_Star, z_MR_Star );
        }

        z_MR_MC.SumScatterFrom( z_MR_Star );
        z = z_MR_MC;
        z.SumScatterUpdate( (T)1, z_MC_Star );
        blas::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        x_MR_Star.FreeAlignments();
        z_MC_Star.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MC,Star> x_MC_Star(g);
        DistMatrix<T,MR,Star> x_MR_Star(g);
        DistMatrix<T,MC,Star> z_MC_Star(g);
        DistMatrix<T,MR,Star> z_MR_Star(g);
        DistMatrix<T,MC,MR  > z(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T,MC,MR  > zTrans(g);

        // Begin the algoritm
        blas::Scal( beta, y );
        x_MC_Star.AlignWith( A );
        x_MR_Star.AlignWith( A );
        z_MC_Star.AlignWith( A );
        z_MR_Star.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_MC_Star.ResizeTo( y.Width(), 1 );
        z_MR_Star.ResizeTo( y.Width(), 1 );
        z_MC_Star.SetToZero();
        z_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        x_MC_Star = x;
        x_MR_Star = x_MC_Star;
        if( shape == Lower )
        {
            blas::internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_Star, x_MR_Star, z_MC_Star, z_MR_Star );
        }
        else
        {
            blas::internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_Star, x_MR_Star, z_MC_Star, z_MR_Star );
        }

        z.SumScatterFrom( z_MC_Star );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( (T)1, z_MR_Star );
        blas::Trans( z_MR_MC, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MC_Star.FreeAlignments();
        x_MR_Star.FreeAlignments();
        z_MC_Star.FreeAlignments();
        z_MR_Star.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(g);
        DistMatrix<T,Star,MR> x_Star_MR(g);
        DistMatrix<T,Star,MC> z_Star_MC(g);
        DistMatrix<T,Star,MR> z_Star_MR(g);
        DistMatrix<T,MC,  MR> z(g);
        DistMatrix<T,MC,  MR> zTrans(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);

        // Begin the algoritm
        blas::Scal( beta, y );
        x_Star_MC.AlignWith( A );
        x_Star_MR.AlignWith( A );
        z_Star_MC.AlignWith( A );
        z_Star_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_Star_MC.ResizeTo( 1, y.Height() );
        z_Star_MR.ResizeTo( 1, y.Height() );
        z_Star_MC.SetToZero();
        z_Star_MR.SetToZero();
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        x_Star_MC = x_Star_MR;
        if( shape == Lower )
        {
            blas::internal::LocalHemvRowAccumulateL
            ( alpha, A, x_Star_MC, x_Star_MR, z_Star_MC, z_Star_MR );
        }
        else
        {
            blas::internal::LocalHemvRowAccumulateU
            ( alpha, A, x_Star_MC, x_Star_MR, z_Star_MC, z_Star_MR );
        }

        z.SumScatterFrom( z_Star_MR );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( (T)1, z_Star_MC );
        blas::Trans( z_MR_MC, zTrans );
        blas::Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        x_Star_MR.FreeAlignments();
        z_Star_MC.FreeAlignments();
        z_Star_MR.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,Star,MC> x_Star_MC(g);
        DistMatrix<T,Star,MR> x_Star_MR(g);
        DistMatrix<T,Star,MC> z_Star_MC(g);
        DistMatrix<T,Star,MR> z_Star_MR(g);
        DistMatrix<T,MC,  MR> z(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);

        // Begin the algoritm
        blas::Scal( beta, y );
        x_Star_MC.AlignWith( A );
        x_Star_MR.AlignWith( A );
        z_Star_MC.AlignWith( A );
        z_Star_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_Star_MC.ResizeTo( 1, y.Width() );
        z_Star_MR.ResizeTo( 1, y.Width() );
        z_Star_MC.SetToZero();
        z_Star_MR.SetToZero();
        //--------------------------------------------------------------------//
        x_Star_MR = x;
        x_Star_MC = x_Star_MR;
        if( shape == Lower )
        {
            blas::internal::LocalHemvRowAccumulateL
            ( alpha, A, x_Star_MC, x_Star_MR, z_Star_MC, z_Star_MR );
        }
        else
        {
            blas::internal::LocalHemvRowAccumulateU
            ( alpha, A, x_Star_MC, x_Star_MR, z_Star_MC, z_Star_MR );
        }

        z_MR_MC.SumScatterFrom( z_Star_MC );
        z = z_MR_MC;
        z.SumScatterUpdate( (T)1, z_Star_MR );
        blas::Axpy( (T)1, z, y );
        //--------------------------------------------------------------------//
        x_Star_MC.FreeAlignments();
        x_Star_MR.FreeAlignments();
        z_Star_MC.FreeAlignments();
        z_Star_MR.FreeAlignments();
        z.FreeAlignments();
        z_MR_MC.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::Hemv
( Shape shape,
  float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& x,
  float beta,        DistMatrix<float,MC,MR>& y );

template void elemental::blas::Hemv
( Shape shape,
  double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& x,
  double beta,        DistMatrix<double,MC,MR>& y );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::Hemv
( Shape shape,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& x,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& y );

template void elemental::blas::Hemv
( Shape shape,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& x,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& y );
#endif

