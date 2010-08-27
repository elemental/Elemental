/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/

namespace {

template<typename T>
void
SetDiagonalToOne( Side side, int offset, DistMatrix<T,MC,MR>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int height = H.Height();
    const int width = H.Width();
    const int localWidth = H.LocalWidth();
    const int r = H.GetGrid().Height();
    const int c = H.GetGrid().Width();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    if( side == Left )
    {
        for( int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const int j = rowShift + jLoc*c;
            const int i = j-offset;     
            if( i >= 0 && i < height && (i-colShift) % r == 0 )
            {
                const int iLoc = (i-colShift)/r;
                H.LocalEntry(iLoc,jLoc) = (T)1;
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
                H.LocalEntry(iLoc,jLoc) = (T)1;
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void 
HalveMainDiagonal( DistMatrix<R,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) /= (R)2;
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
FixDiagonal
( const DistMatrix<complex<R>,Star,Star>& t,
        DistMatrix<complex<R>,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("FixDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) = complex<R>(1)/t.LocalEntry(j,0);
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

}

