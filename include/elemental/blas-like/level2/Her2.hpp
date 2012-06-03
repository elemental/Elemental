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

template<typename T>
inline void
Her2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
           const DistMatrix<T>& y,
                 DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her2");
    if( A.Grid() != x.Grid() || x.Grid() != y.Grid() )
        throw std::logic_error
        ("{A,x,y} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        std::ostringstream msg;
        msg << "A must conform with x: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n"
            << "  y ~ " << y.Height() << " x " << y.Width() << "\n";
        throw std::logic_error( msg.str() );
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
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,MC,STAR> y_MC_STAR(g);
        DistMatrix<T,MR,STAR> y_MR_STAR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        y_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        y_MC_STAR = y;
        y_MR_STAR = y_MC_STAR;

        if( uplo == LOWER )
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
                      (                 x_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(y_MR_STAR.GetLocalEntry(jLoc,0)) + 
                             y_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(x_MR_STAR.GetLocalEntry(jLoc,0)) ) );
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
                      (                 x_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(y_MR_STAR.GetLocalEntry(jLoc,0)) + 
                             y_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(x_MR_STAR.GetLocalEntry(jLoc,0)) ) );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
        y_MC_STAR.FreeAlignments();
        y_MR_STAR.FreeAlignments();
    }
    else if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);
        DistMatrix<T,STAR,MC> y_STAR_MC(g);
        DistMatrix<T,STAR,MR> y_STAR_MR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        y_STAR_MC.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;
        y_STAR_MR = y;
        y_STAR_MC = y_STAR_MR;

        if( uplo == LOWER )
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
                      (                 x_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(y_STAR_MR.GetLocalEntry(0,jLoc)) + 
                             y_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(x_MR_STAR.GetLocalEntry(jLoc,0)) ) );
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
                      (                 x_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(y_STAR_MR.GetLocalEntry(0,jLoc)) +
                             y_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(x_MR_STAR.GetLocalEntry(jLoc,0)) ) );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
        y_STAR_MC.FreeAlignments();
        y_STAR_MR.FreeAlignments();
    }
    else if( y.Width() == 1 )
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        DistMatrix<T,MC,STAR> y_MC_STAR(g);
        DistMatrix<T,MR,STAR> y_MR_STAR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        y_MC_STAR.AlignWith( A );
        y_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        y_MC_STAR = y;
        y_MR_STAR = y_MC_STAR;

        if( uplo == LOWER )
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
                      (                 x_STAR_MC.GetLocalEntry(0,iLoc)* 
                        Conj(y_MR_STAR.GetLocalEntry(jLoc,0)) + 
                             y_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(x_STAR_MR.GetLocalEntry(0,jLoc)) ) );
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
                      (                 x_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(y_MR_STAR.GetLocalEntry(jLoc,0)) +
                             y_MC_STAR.GetLocalEntry(iLoc,0)*
                        Conj(x_STAR_MR.GetLocalEntry(0,jLoc)) ) );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
        y_MC_STAR.FreeAlignments();
        y_MR_STAR.FreeAlignments();
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);
        DistMatrix<T,STAR,MC> y_STAR_MC(g);
        DistMatrix<T,STAR,MR> y_STAR_MR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        y_STAR_MC.AlignWith( A );
        y_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;
        y_STAR_MR = y;
        y_STAR_MC = y_STAR_MR;

        if( uplo == LOWER )
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
                      (                 x_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(y_STAR_MR.GetLocalEntry(0,jLoc)) + 
                             y_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(x_STAR_MR.GetLocalEntry(0,jLoc)) ) );
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
                      (                 x_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(y_STAR_MR.GetLocalEntry(0,jLoc)) + 
                             y_STAR_MC.GetLocalEntry(0,iLoc)*
                        Conj(x_STAR_MR.GetLocalEntry(0,jLoc)) ) );
                }
            }
        }
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
        y_STAR_MC.FreeAlignments();
        y_STAR_MR.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
