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
void
SetDiagonalToOne( LeftOrRight side, int offset, DistMatrix<T,MC,MR>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int height = H.Height();
    const int width = H.Width();
    const int localWidth = H.LocalWidth();
    const int r = H.Grid().Height();
    const int c = H.Grid().Width();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    if( side == LEFT )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            const int i = j-offset;     
            if( i >= 0 && i < height && (i-colShift) % r == 0 )
            {
                const int iLoc = (i-colShift)/r;
                H.SetLocalEntry(iLoc,jLoc,1);
            }
        }
    }
    else
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            const int i = j-offset+height-width;
            if( i >= 0 && i < height && (i-colShift) % r == 0 )
            {
                const int iLoc = (i-colShift)/r;
                H.SetLocalEntry(iLoc,jLoc,1);
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void 
HalveMainDiagonal( DistMatrix<R,STAR,STAR>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
    {
        const R value = SInv.GetLocalEntry(j,j);
        SInv.SetLocalEntry(j,j,value/2);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void
FixDiagonal
( Conjugation conjugation,
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R>,STAR,STAR>& SInv )
{
#ifndef RELEASE
    PushCallStack("FixDiagonal");
#endif
    if( conjugation == CONJUGATED )
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const Complex<R> value = Complex<R>(1)/t.GetLocalEntry(j,0);
            SInv.SetLocalEntry(j,j,value);
        }
    }
    else
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const Complex<R> value = Complex<R>(1)/Conj(t.GetLocalEntry(j,0));
            SInv.SetLocalEntry(j,j,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // internal
} // elem
