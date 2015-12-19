/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LATTICE_HPP
#define EL_LATTICE_HPP

namespace El {

// TODO: Maintain B in BigInt form

struct LLLInfo
{
    Int nullity; 
    Int numSwaps;
};

template<typename Real>
struct LLLCtrl
{
    Real delta=Real(3)/Real(4);

    // A 'weak' LLL reduction only ensures that | R(i,i+1) / R(i,i) | is
    // bounded above by one-half (for complex data, by sqrt(2)/2)
    bool weak=false;

    // Preprocessing with a "rank-obscuring" column-pivoted QR factorization
    // (in the manner suggested by Wubben et al.) tends to greatly decrease
    // the number of swaps within LLL
    bool presort=true;
    bool smallestFirst=true;

    // If the size-reduced column has a squared two-norm that is less than or
    // equal to `reorthogTol` times the square of its original two-norm, then
    // it is reorthogonalized
    Real reorthogTol=0;

    // If a size-reduced column has a two-norm less than or equal to 'zeroTol',
    // then it is interpreted as a zero vector (and forced to zero)
    Real zeroTol=limits::Epsilon<Real>();

    bool progress=false;
    bool time=false;
};

template<typename F>
LLLInfo LLL
( Matrix<F>& B,
  Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
Base<F> LLLDelta
( const Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
void LatticeGramSchmidt( const Matrix<F>& B, Matrix<F>& G, Matrix<F>& M );

} // namespace El

#endif // ifndef EL_LATTICE_HPP
