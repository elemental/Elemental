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

template<typename F>
Int LLL
( Matrix<F>& B,
  Matrix<F>& QR,
  Base<F> delta,
  Base<F> innerTol=0,
  bool presort=false,
  bool smallestFirst=true,
  bool progress=false,
  bool time=false );
template<typename F>
Base<F> LLLDelta( const Matrix<F>& QR );

template<typename F>
void LatticeGramSchmidt( const Matrix<F>& B, Matrix<F>& G, Matrix<F>& M );

} // namespace El

#endif // ifndef EL_LATTICE_HPP
