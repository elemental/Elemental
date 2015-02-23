/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EQUILIBRATE_HPP
#define EL_EQUILIBRATE_HPP

namespace El {

// Geometric rescaling (ala Fourer, which led to Saunders's gmscale.m)
// ===================================================================
template<typename F>
void GeomEquil
( Matrix<F>& A, 
  Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void GeomEquil
( AbstractDistMatrix<F>& A, 
  AbstractDistMatrix<Base<F>>& dRow, AbstractDistMatrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void GeomEquil
( SparseMatrix<F>& A, 
  Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void GeomEquil
( DistSparseMatrix<F>& A, 
  DistMultiVec<Base<F>>& dRow, DistMultiVec<Base<F>>& dCol,
  bool progress=false );

} // namespace El

#endif // ifndef EL_EQUILIBRATE_HPP
