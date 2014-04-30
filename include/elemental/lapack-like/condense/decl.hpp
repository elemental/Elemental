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

namespace HermitianTridiagApproachNS {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace HermitianTridiagApproachNS;

struct HermitianTridiagCtrl {
    HermitianTridiagApproach approach;
    GridOrder order;

    HermitianTridiagCtrl() 
    : approach(HERMITIAN_TRIDIAG_SQUARE), order(ROW_MAJOR) 
    { }
};

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t );
template<typename F>
void HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t,
  const HermitianTridiagCtrl ctrl=HermitianTridiagCtrl() );

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A,
  const HermitianTridiagCtrl ctrl=HermitianTridiagCtrl() );

} // namespace elem

#endif // ifndef ELEM_CONDENSE_DECL_HPP
