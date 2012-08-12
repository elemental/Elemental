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
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
        throw std::logic_error("x must conform with A");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    blas::Her
    ( uploChar, m, alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
                 DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
    if( A.Grid() != x.Grid() )
        throw std::logic_error("{A,x} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( A.Height() != xLength )
    {
        std::ostringstream msg;
        msg << "A must conform with x: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x ~ " << x.Height() << " x " << x.Width() << "\n";
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

    if( x.Width() == 1 )
    {
        DistMatrix<T,MC,STAR> x_MC_STAR(g);
        DistMatrix<T,MR,STAR> x_MR_STAR(g);

        x_MC_STAR.AlignWith( A );
        x_MR_STAR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_MC_STAR = x;
        x_MR_STAR = x_MC_STAR;

        const T* xLocal = x_MC_STAR.LockedLocalBuffer();
        if( uplo == LOWER )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);

                const T gamma = alpha*Conj(x_MR_STAR.GetLocal(jLocal,0));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=heightAboveDiag; iLocal<localHeight; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal];
            }
        }
        else
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);

                const T gamma = alpha*Conj(x_MR_STAR.GetLocal(jLocal,0));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=0; iLocal<heightToDiag; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal];
            }
        }
        //--------------------------------------------------------------------//
        x_MC_STAR.FreeAlignments();
        x_MR_STAR.FreeAlignments();
    }
    else
    {
        DistMatrix<T,STAR,MC> x_STAR_MC(g);
        DistMatrix<T,STAR,MR> x_STAR_MR(g);

        x_STAR_MC.AlignWith( A );
        x_STAR_MR.AlignWith( A );
        //--------------------------------------------------------------------//
        x_STAR_MR = x;
        x_STAR_MC = x_STAR_MR;

        const T* xLocal = x_STAR_MC.LockedLocalBuffer();
        const int incx = x_STAR_MC.LocalLDim();
        if( uplo == LOWER )
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightAboveDiag = LocalLength(j,colShift,r);

                const T gamma = alpha*Conj(x_STAR_MR.GetLocal(0,jLocal));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=heightAboveDiag; iLocal<localHeight; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal*incx];
            }
        }
        else
        {
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
            {
                const int j = rowShift + jLocal*c;
                const int heightToDiag = LocalLength(j+1,colShift,r);

                const T gamma = alpha*Conj(x_STAR_MR.GetLocal(0,jLocal));
                T* ALocalCol = A.LocalBuffer(0,jLocal);
                for( int iLocal=0; iLocal<heightToDiag; ++iLocal )
                    ALocalCol[iLocal] += gamma*xLocal[iLocal*incx];
            }
        }
        //--------------------------------------------------------------------//
        x_STAR_MC.FreeAlignments();
        x_STAR_MR.FreeAlignments();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
