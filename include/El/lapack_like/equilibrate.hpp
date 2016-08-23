/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_EQUILIBRATE_HPP
#define EL_EQUILIBRATE_HPP

namespace El {

// Ruiz scaling
// ============
template<typename F>
void RuizEquil
( Matrix<F>& A,
  Matrix<Base<F>>& dRow,
  Matrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void RuizEquil
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& dRow,
  ElementalMatrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void RuizEquil
( SparseMatrix<F>& A,
  Matrix<Base<F>>& dRow,
  Matrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void RuizEquil
( DistSparseMatrix<F>& A,
  DistMultiVec<Base<F>>& dRow,
  DistMultiVec<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void StackedRuizEquil
( Matrix<F>& A,
  Matrix<F>& B,
  Matrix<Base<F>>& dRowA,
  Matrix<Base<F>>& dRowB,
  Matrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void StackedRuizEquil
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& B,
  ElementalMatrix<Base<F>>& dRowA,
  ElementalMatrix<Base<F>>& dRowB,
  ElementalMatrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void StackedRuizEquil
( SparseMatrix<F>& A,
  SparseMatrix<F>& B,
  Matrix<Base<F>>& dRowA,
  Matrix<Base<F>>& dRowB,
  Matrix<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void StackedRuizEquil
( DistSparseMatrix<F>& A,
  DistSparseMatrix<F>& B,
  DistMultiVec<Base<F>>& dRowA,
  DistMultiVec<Base<F>>& dRowB,
  DistMultiVec<Base<F>>& dCol,
  bool progress=false );

template<typename F>
void SymmetricRuizEquil
( Matrix<F>& A,
  Matrix<Base<F>>& d,
  Int maxiter=3, bool progress=false );

template<typename F>
void SymmetricRuizEquil
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& d,
  Int maxiter=3, bool progress=false );

template<typename F>
void SymmetricRuizEquil
( SparseMatrix<F>& A,
  Matrix<Base<F>>& d,
  Int maxiter=3, bool progress=false );

template<typename F>
void SymmetricRuizEquil
( DistSparseMatrix<F>& A,
  DistMultiVec<Base<F>>& d,
  Int maxiter=3, bool progress=false );

// Geometric rescaling (ala Fourer, which led to Saunders's gmscale.m)
// ===================================================================
template<typename F>
void GeomEquil
( Matrix<F>& A, 
  Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void GeomEquil
( ElementalMatrix<F>& A, 
  ElementalMatrix<Base<F>>& dRow, ElementalMatrix<Base<F>>& dCol,
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

// NOTE: It is assumed that the matrix is explicitly symmetric/Hermitian
//       (either will be equivalent for this routine)
template<typename F>
void SymmetricGeomEquil
( Matrix<F>& A, Matrix<Base<F>>& d, bool progress=false );
template<typename F>
void SymmetricGeomEquil
( ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& d, 
  bool progress=false );
template<typename F>
void SymmetricGeomEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& d, bool progress=false );
template<typename F>
void SymmetricGeomEquil
( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, bool progress=false );

template<typename F>
void StackedGeomEquil
( Matrix<F>& A, Matrix<F>& B,
  Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB, Matrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void StackedGeomEquil
( ElementalMatrix<F>& A, ElementalMatrix<F>& B,
  ElementalMatrix<Base<F>>& dRowA, 
  ElementalMatrix<Base<F>>& dRowB, 
  ElementalMatrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void StackedGeomEquil
( SparseMatrix<F>& A, SparseMatrix<F>& B,
  Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB, 
  Matrix<Base<F>>& dCol,
  bool progress=false );
template<typename F>
void StackedGeomEquil
( DistSparseMatrix<F>& A, DistSparseMatrix<F>& B,
  DistMultiVec<Base<F>>& dRowA, DistMultiVec<Base<F>>& dRowB, 
  DistMultiVec<Base<F>>& dCol,
  bool progress=false );

// Diagonal equilibration
// ======================
template<typename F>
void SymmetricDiagonalEquil
( Matrix<F>& A, Matrix<Base<F>>& d, bool progress=false );
template<typename F>
void SymmetricDiagonalEquil
( ElementalMatrix<F>& A, ElementalMatrix<Base<F>>& d, 
  bool progress=false );
template<typename F>
void SymmetricDiagonalEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& d, bool progress=false );
template<typename F>
void SymmetricDiagonalEquil
( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, 
  bool progress=false, bool time=false );

} // namespace El

#endif // ifndef EL_EQUILIBRATE_HPP
