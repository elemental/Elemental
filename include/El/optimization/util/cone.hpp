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
template<typename Field>
void Broadcast
(       Matrix<Field>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Field>
void Broadcast
(       AbstractDistMatrix<Field>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Field>
void Broadcast
(       DistMultiVec<Field>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// AllReduce
// =========
// Fill each subcone with the reduction over each cone
template<typename Field>
void AllReduce
(       Matrix<Field>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  mpi::Op op=mpi::SUM );
template<typename Field>
void AllReduce
(       AbstractDistMatrix<Field>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  mpi::Op op=mpi::SUM, Int cutoff=1000 );
template<typename Field>
void AllReduce
(       DistMultiVec<Field>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  mpi::Op op=mpi::SUM, Int cutoff=1000 );

// A specialization of Ruiz scaling which respects a product of cones
// ==================================================================
template<typename Field>
void RuizEquil
(       Matrix<Field>& A,
        Matrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename Field>
void RuizEquil
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Base<Field>>& dRowA,
        AbstractDistMatrix<Base<Field>>& dRowB,
        AbstractDistMatrix<Base<Field>>& dCol,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename Field>
void RuizEquil
(       SparseMatrix<Field>& A,
        SparseMatrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename Field>
void RuizEquil
(       DistSparseMatrix<Field>& A,
        DistSparseMatrix<Field>& B,
        DistMultiVec<Base<Field>>& dRowA,
        DistMultiVec<Base<Field>>& dRowB,
        DistMultiVec<Base<Field>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

// A specialization of GeomEquil which respects a product of cones
// ===============================================================
template<typename Field>
void GeomEquil
(       Matrix<Field>& A,
        Matrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename Field>
void GeomEquil
(       AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Field>& BPre,
        AbstractDistMatrix<Base<Field>>& dRowAPre,
        AbstractDistMatrix<Base<Field>>& dRowBPre,
        AbstractDistMatrix<Base<Field>>& dColPre,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename Field>
void GeomEquil
(       SparseMatrix<Field>& A,
        SparseMatrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename Field>
void GeomEquil
(       DistSparseMatrix<Field>& A,
        DistSparseMatrix<Field>& B,
        DistMultiVec<Base<Field>>& dRowA,
        DistMultiVec<Base<Field>>& dRowB,
        DistMultiVec<Base<Field>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

} // namespace cone
} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_CONE_HPP
