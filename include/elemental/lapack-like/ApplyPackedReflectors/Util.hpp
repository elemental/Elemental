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

template<typename F> 
inline void
FixDiagonal
( Conjugation conjugation, const Matrix<F>& t, Matrix<F>& SInv )
{
#ifndef RELEASE
    CallStackEntry cse("FixDiagonal");
#endif
    if( conjugation == CONJUGATED )
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const F value = F(1)/Conj(t.Get(j,0));
            SInv.Set(j,j,value);
        }
    }
    else
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const F value = F(1)/t.Get(j,0);
            SInv.Set(j,j,value);
        }
    }
}

template<typename F> 
inline void
FixDiagonal
( Conjugation conjugation,
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F,STAR,STAR>& SInv )
{
#ifndef RELEASE
    CallStackEntry cse("FixDiagonal");
#endif
    if( conjugation == CONJUGATED )
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const F value = F(1)/Conj(t.GetLocal(j,0));
            SInv.SetLocal(j,j,value);
        }
    }
    else
    {
        for( int j=0; j<SInv.Height(); ++j )
        {
            const F value = F(1)/t.GetLocal(j,0);
            SInv.SetLocal(j,j,value);
        }
    }
}

} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_UTIL_HPP
