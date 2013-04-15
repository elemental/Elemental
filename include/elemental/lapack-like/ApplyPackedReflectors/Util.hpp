/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_UTIL_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_UTIL_HPP

namespace elem {

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

} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_UTIL_HPP
