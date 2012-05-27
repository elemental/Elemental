/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename T>
inline void
GemvN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y )
{
#ifndef RELEASE
    PushCallStack("internal::GemvN");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( ( x.Width() != 1 && x.Height() != 1 ) ||
        ( y.Width() != 1 && y.Height() != 1 )   )
        throw std::logic_error("x and y are assumed to be vectors");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != yLength || A.Width() != xLength )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemvN: \n"
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
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);

        // Start the algorithm
        Scal( beta, y );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_MR_STAR = x;
        Gemv
        ( NORMAL,
          alpha, A.LockedLocalMatrix(), 
                 x_MR_STAR.LockedLocalMatrix(),
          (T)0,  z_MC_STAR.LocalMatrix() );
        y.SumScatterUpdate( (T)1, z_MC_STAR );
        //--------------------------------------------------------------------//
        x_MR_STAR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> z_MC_STAR(g);
        DistMatrix<T,MC,MR  > z(g);
        DistMatrix<T,MC,MR  > zTrans(g);

        // Start the algorithm
        Scal( beta, y );
        x_MR_STAR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_MR_STAR = x;
        Gemv
        ( NORMAL,
          alpha, A.LockedLocalMatrix(),
                 x_MR_STAR.LockedLocalMatrix(),
          (T)0,  z_MC_STAR.LocalMatrix() );
        z.SumScatterFrom( z_MC_STAR );
        Transpose( z, zTrans );
        Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_MR_STAR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
        z.FreeAlignments();
        zTrans.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        // Temporary distributions
        DistMatrix<T,STAR,MR  > x_STAR_MR(g);
        DistMatrix<T,MC,  STAR> z_MC_STAR(g);

        // Start the algorithm
        Scal( beta, y );
        x_STAR_MR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.ResizeTo( A.Height(), 1 );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        Gemv
        ( NORMAL,
          alpha, A.LockedLocalMatrix(), 
                 x_STAR_MR.LockedLocalMatrix(),
          (T)0,  z_MC_STAR.LocalMatrix() );
        y.SumScatterUpdate( (T)1, z_MC_STAR );
        //--------------------------------------------------------------------//
        x_STAR_MR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
    }
    else
    {
        // Temporary distributions
        DistMatrix<T,STAR,MR  > x_STAR_MR(g);
        DistMatrix<T,MC,  STAR> z_MC_STAR(g);
        DistMatrix<T,MC,  MR  > z(g);
        DistMatrix<T,MC,  MR  > zTrans(g);

        // Start the algorithm
        Scal( beta, y );
        x_STAR_MR.AlignWith( A );
        z_MC_STAR.AlignWith( A );
        z_MC_STAR.ResizeTo( A.Height(), 1 );
        z.AlignWith( y );
        zTrans.AlignWith( y );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        Gemv
        ( NORMAL,
          alpha, A.LockedLocalMatrix(),
                 x_STAR_MR.LockedLocalMatrix(),
          (T)0,  z_MC_STAR.LocalMatrix() );
        z.SumScatterFrom( z_MC_STAR );
        Transpose( z, zTrans );
        Axpy( (T)1, zTrans, y );
        //--------------------------------------------------------------------//
        x_STAR_MR.FreeAlignments();
        z_MC_STAR.FreeAlignments();
        z.FreeAlignments();
        zTrans.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
