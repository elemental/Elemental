/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_UTIL_CONE_HPP
#define EL_OPTIMIZATION_UTIL_CONE_HPP

namespace El {
namespace cone {

// Broadcast
// =========
// Replicate the entry in the root position in each cone over the entire cone
template<typename F>
void Broadcast
(       Matrix<F>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename F>
void Broadcast
(       ElementalMatrix<F>& x, 
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename F>
void Broadcast
(       DistMultiVec<F>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// AllReduce
// =========
// Fill each subcone with the reduction over each cone
template<typename F>
void AllReduce
(       Matrix<F>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  mpi::Op op=mpi::SUM );
template<typename F>
void AllReduce
(       ElementalMatrix<F>& x, 
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds, 
  mpi::Op op=mpi::SUM, Int cutoff=1000 );
template<typename F>
void AllReduce
(       DistMultiVec<F>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, 
  mpi::Op op=mpi::SUM, Int cutoff=1000 );

// A specialization of Ruiz scaling which respects a product of cones
// ==================================================================
template<typename F>
void RuizEquil
(       Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void RuizEquil
(       ElementalMatrix<F>& A,
        ElementalMatrix<F>& B,
        ElementalMatrix<Base<F>>& dRowA,
        ElementalMatrix<Base<F>>& dRowB,
        ElementalMatrix<Base<F>>& dCol,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename F>
void RuizEquil
(       SparseMatrix<F>& A,
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void RuizEquil
(       DistSparseMatrix<F>& A,
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA,
        DistMultiVec<Base<F>>& dRowB,
        DistMultiVec<Base<F>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

// A specialization of GeomEquil which respects a product of cones
// ===============================================================
template<typename F>
void GeomEquil
(       Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void GeomEquil
(       ElementalMatrix<F>& APre,
        ElementalMatrix<F>& BPre,
        ElementalMatrix<Base<F>>& dRowAPre,
        ElementalMatrix<Base<F>>& dRowBPre,
        ElementalMatrix<Base<F>>& dColPre,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename F>
void GeomEquil
(       SparseMatrix<F>& A,
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void GeomEquil
(       DistSparseMatrix<F>& A,
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA,
        DistMultiVec<Base<F>>& dRowB,
        DistMultiVec<Base<F>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

} // namespace cone
} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_CONE_HPP
