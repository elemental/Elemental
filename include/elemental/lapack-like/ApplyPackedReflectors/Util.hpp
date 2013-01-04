/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T> 
inline void
SetDiagonalToOne( LeftOrRight side, int offset, Matrix<T>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int height = H.Height();
    const int width = H.Width();

    if( side == LEFT )
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset;     
            if( i >= 0 && i < height )
                H.Set(i,j,1);
        }
    }
    else
    {
        for( int j=0; j<width; ++j )
        {
            const int i = j-offset+height-width;
            if( i >= 0 && i < height )
                H.Set(i,j,1);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> 
inline void
SetDiagonalToOne( LeftOrRight side, int offset, DistMatrix<T>& H )
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
                H.SetLocal(iLoc,jLoc,1);
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
                H.SetLocal(iLoc,jLoc,1);
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void 
HalveMainDiagonal( Matrix<R>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
    {
        const R value = SInv.Get(j,j);
        SInv.Set(j,j,value/2);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void 
HalveMainDiagonal( DistMatrix<R,STAR,STAR>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
    {
        const R value = SInv.GetLocal(j,j);
        SInv.SetLocal(j,j,value/2);
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
FixDiagonal
( Conjugation conjugation,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& SInv )
{
#ifndef RELEASE
    PushCallStack("FixDiagonal");
#endif
    if( conjugation == CONJUGATED )
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const Complex<R> value = Complex<R>(1)/Conj(t.Get(j,0));
            SInv.Set(j,j,value);
        }
    }
    else
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const Complex<R> value = Complex<R>(1)/t.Get(j,0);
            SInv.Set(j,j,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
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
            const Complex<R> value = Complex<R>(1)/Conj(t.GetLocal(j,0));
            SInv.SetLocal(j,j,value);
        }
    }
    else
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const Complex<R> value = Complex<R>(1)/t.GetLocal(j,0);
            SInv.SetLocal(j,j,value);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // internal
} // elem
