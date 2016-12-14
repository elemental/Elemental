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
template<typename Field>
void RuizEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void RuizEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& dRow,
  AbstractDistMatrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void RuizEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void RuizEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& dRow,
  DistMultiVec<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void StackedRuizEquil
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void StackedRuizEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Base<Field>>& dRowA,
  AbstractDistMatrix<Base<Field>>& dRowB,
  AbstractDistMatrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void StackedRuizEquil
( SparseMatrix<Field>& A,
  SparseMatrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void StackedRuizEquil
( DistSparseMatrix<Field>& A,
  DistSparseMatrix<Field>& B,
  DistMultiVec<Base<Field>>& dRowA,
  DistMultiVec<Base<Field>>& dRowB,
  DistMultiVec<Base<Field>>& dCol,
  bool progress=false );

template<typename Field>
void SymmetricRuizEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& d,
  Int maxiter=3, bool progress=false );

template<typename Field>
void SymmetricRuizEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& d,
  Int maxiter=3, bool progress=false );

template<typename Field>
void SymmetricRuizEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& d,
  Int maxiter=3, bool progress=false );

template<typename Field>
void SymmetricRuizEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& d,
  Int maxiter=3, bool progress=false );

// Geometric rescaling (ala Fourer, which led to Saunders's gmscale.m)
// ===================================================================
template<typename Field>
void GeomEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void GeomEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& dRow,
  AbstractDistMatrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void GeomEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void GeomEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& dRow,
  DistMultiVec<Base<Field>>& dCol,
  bool progress=false );

// NOTE: It is assumed that the matrix is explicitly symmetric/Hermitian
//       (either will be equivalent for this routine)
template<typename Field>
void SymmetricGeomEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricGeomEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricGeomEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricGeomEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& d,
  bool progress=false );

template<typename Field>
void StackedGeomEquil
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void StackedGeomEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Base<Field>>& dRowA,
  AbstractDistMatrix<Base<Field>>& dRowB,
  AbstractDistMatrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void StackedGeomEquil
( SparseMatrix<Field>& A,
  SparseMatrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
  bool progress=false );
template<typename Field>
void StackedGeomEquil
( DistSparseMatrix<Field>& A,
  DistSparseMatrix<Field>& B,
  DistMultiVec<Base<Field>>& dRowA,
  DistMultiVec<Base<Field>>& dRowB,
  DistMultiVec<Base<Field>>& dCol,
  bool progress=false );

// Diagonal equilibration
// ======================
template<typename Field>
void SymmetricDiagonalEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricDiagonalEquil
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricDiagonalEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& d,
  bool progress=false );
template<typename Field>
void SymmetricDiagonalEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& d,
  bool progress=false, bool time=false );

} // namespace El

#endif // ifndef EL_EQUILIBRATE_HPP
