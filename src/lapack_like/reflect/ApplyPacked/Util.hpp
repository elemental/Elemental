/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_UTIL_HPP
#define EL_APPLYPACKEDREFLECTORS_UTIL_HPP

namespace El {

template<typename F> 
void FixDiagonal
( Conjugation conjugation, const Matrix<F>& t, Matrix<F>& SInv )
{
    DEBUG_CSE
    for( Int j=0; j<SInv.Height(); ++j )
    {
        const F value = t(j);
        if( conjugation == CONJUGATED )
            SInv(j,j) = F(1)/Conj(value);
        else
            SInv(j,j) = F(1)/value;
    }
}

template<typename F> 
void FixDiagonal
( Conjugation conjugation,
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F,STAR,STAR>& SInv )
{
    DEBUG_CSE
    auto& tLoc = t.LockedMatrix();
    auto& SInvLoc = SInv.Matrix();
    for( Int j=0; j<SInv.Height(); ++j )
    {
        const F value = tLoc(j);
        if( conjugation == CONJUGATED )
            SInvLoc(j,j) = F(1)/Conj(value);
        else
            SInvLoc(j,j) = F(1)/value;
    }
}

} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_UTIL_HPP
