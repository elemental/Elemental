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

#include "./Hemv/HemvL.hpp"
#include "./Hemv/HemvU.hpp"

template<typename T>
inline void
elemental::basic::Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("basic::Hemv");
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
        DistMatrix<T,MC,MR  > z(g);

        // Begin the algoritm
        basic::Scal( beta, y );
        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z.AlignWith( y );
        z_MC_STAR.ResizeTo( y.Height(), 1 );
        z_MR_STAR.ResizeTo( y.Height(), 1 );
        z_MC_STAR.SetToZero();
        z_MR_STAR.SetToZero();
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            basic::internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }
        else
        {
            basic::internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }

        z_MR_MC.SumScatterFrom( z_MR_STAR );
        z = z_MR_MC;
        z.SumScatterUpdate( (T)1, z_MC_STAR );
        basic::Axpy( (T)1, z, y );
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
        DistMatrix<T,MC,MR  > z(g);
        DistMatrix<T,MR,MC  > z_MR_MC(g);
        DistMatrix<T,MC,MR  > zTrans(g);

        // Begin the algoritm
        basic::Scal( beta, y );
        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MR_STAR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_MC_STAR.ResizeTo( y.Width(), 1 );
        z_MR_STAR.ResizeTo( y.Width(), 1 );
        z_MC_STAR.SetToZero();
        z_MR_STAR.SetToZero();
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        if( uplo == LOWER )
        {
            basic::internal::LocalHemvColAccumulateL
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }
        else
        {
            basic::internal::LocalHemvColAccumulateU
            ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        }

        z.SumScatterFrom( z_MC_STAR );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( (T)1, z_MR_STAR );
        basic::Transpose( z_MR_MC, zTrans );
        basic::Axpy( (T)1, zTrans, y );
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
        DistMatrix<T,MC,  MR> z(g);
        DistMatrix<T,MC,  MR> zTrans(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);

        // Begin the algoritm
        basic::Scal( beta, y );
        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_STAR_MC.ResizeTo( 1, y.Height() );
        z_STAR_MR.ResizeTo( 1, y.Height() );
        z_STAR_MC.SetToZero();
        z_STAR_MR.SetToZero();
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            basic::internal::LocalHemvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }
        else
        {
            basic::internal::LocalHemvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }

        z.SumScatterFrom( z_STAR_MR );
        z_MR_MC = z;
        z_MR_MC.SumScatterUpdate( (T)1, z_STAR_MC );
        basic::Transpose( z_MR_MC, zTrans );
        basic::Axpy( (T)1, zTrans, y );
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
        DistMatrix<T,MC,  MR> z(g);
        DistMatrix<T,MR,  MC> z_MR_MC(g);

        // Begin the algoritm
        basic::Scal( beta, y );
        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        z_STAR_MC.AlignWith( A );
        z_STAR_MR.AlignWith( A );
        z.AlignWith( y );
        z_MR_MC.AlignWith( y );
        z_STAR_MC.ResizeTo( 1, y.Width() );
        z_STAR_MR.ResizeTo( 1, y.Width() );
        z_STAR_MC.SetToZero();
        z_STAR_MR.SetToZero();
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        if( uplo == LOWER )
        {
            basic::internal::LocalHemvRowAccumulateL
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }
        else
        {
            basic::internal::LocalHemvRowAccumulateU
            ( alpha, A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        }

        z_MR_MC.SumScatterFrom( z_STAR_MC );
        z = z_MR_MC;
        z.SumScatterUpdate( (T)1, z_STAR_MR );
        basic::Axpy( (T)1, z, y );
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
