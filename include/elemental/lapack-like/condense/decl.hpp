/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONDENSE_DECL_HPP
#define ELEM_CONDENSE_DECL_HPP

namespace elem {

namespace hermitian_tridiag_approach_wrapper {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace hermitian_tridiag_approach_wrapper;

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t );
template<typename F>
void HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t );

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HermitianTridiag( UpperOrLower uplo, DistMatrix<F>& A );

void SetHermitianTridiagApproach( HermitianTridiagApproach approach );
HermitianTridiagApproach GetHermitianTridiagApproach();

// If dropping down to a square grid, the two simplest approaches are to take 
// the first r^2 processes from the original grid (for an r x r grid) and to
// either order them column-major or row-major to form the square grid.
void SetHermitianTridiagGridOrder( GridOrder order );
GridOrder GetHermitianTridiagGridOrder();

} // namespace elem

#endif // ifndef ELEM_CONDENSE_DECL_HPP
