/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_INVERSE_HERMITIAN_HPP
#define ELEM_INVERSE_HERMITIAN_HPP

#include "./Symmetric.hpp"

namespace elem {

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
inline void
LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHermitianInverse"))
    SymmetricInverse( uplo, A.Matrix(), true, pivotType );
}

} // namespace elem

#endif // ifndef ELEM_INVERSE_HERMITIAN_HPP
