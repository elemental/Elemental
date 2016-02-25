/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PERM_HPP
#define EL_PERM_HPP

#include "perm/Permutation.hpp"
#include "perm/DistPermutation.hpp"

namespace El {

// Convert a pivot sequence to a partial permutation vector
// ========================================================
// NOTE: These routine are now deprecated
void PivotsToPartialPermutation
( const Matrix<Int>& pivots,
        Matrix<Int>& p,
        Matrix<Int>& pInv,
        Int offset=0 );
void PivotsToPartialPermutation
( const AbstractDistMatrix<Int>& pivots,
        AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv,
        Int offset=0 );
void PivotsToPartialPermutation
( const DistMatrix<Int,STAR,STAR>& pivots,
        AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv,
        Int offset=0 );

} // namespace El

#endif // ifndef EL_PERM_HPP
