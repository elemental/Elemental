/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
inline Real DampScaling( Real alpha )
{
    const Real tol = Pow(Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else 
        return Max(alpha,tol);
}

template<typename F>
void ConeRuizEquil
(       Matrix<F>& A, 
        Matrix<F>& B, 
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB, 
        Matrix<Base<F>>& dCol, 
  const Matrix<Int>& orders,  
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeRuizEquil"))
    LogicError("This routine is not yet written");
}

template<typename F>
void ConeRuizEquil
(       AbstractDistMatrix<F>& APre, 
        AbstractDistMatrix<F>& BPre,
        AbstractDistMatrix<Base<F>>& dRowAPre, 
        AbstractDistMatrix<Base<F>>& dRowBPre,
        AbstractDistMatrix<Base<F>>& dColPre,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeRuizEquil"))
    typedef Base<F> Real;

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr     = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto BPtr     = ReadWriteProxy<F,MC,MR>(&BPre,control);
    auto dRowAPtr = WriteProxy<Real,MC,STAR>(&dRowAPre,control); 
    auto dRowBPtr = WriteProxy<Real,MC,STAR>(&dRowBPre,control); 
    auto dColPtr  = WriteProxy<Real,MR,STAR>(&dColPre,control);
    auto& A = *APtr;
    auto& B = *BPtr;
    auto& dRowA = *dRowAPtr;
    auto& dRowB = *dRowBPtr;
    auto& dCol  = *dColPtr;

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    const Int nLocal = A.LocalWidth();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4; 

    DistMatrix<Real,MC,STAR> rowScale(A.Grid());
    DistMatrix<Real,MR,STAR> colScale(A.Grid()), colScaleB(B.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        ColumnMaxNorms( B, colScaleB );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            colScale.SetLocal
            ( jLoc, 0, Max(colScale.GetLocal(jLoc,0),
                           colScaleB.GetLocal(jLoc,0)) );
        EntrywiseMap( colScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, rowScale );
        EntrywiseMap( rowScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScale, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScale, A );

        RowMaxNorms( B, rowScale );
        ConeAllReduce( rowScale, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( rowScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScale, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScale, B );
    }
    SetIndent( indent );
}

template<typename F>
void ConeRuizEquil
(       SparseMatrix<F>& A, 
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeRuizEquil"))
    typedef Base<F> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    // For now, simply hard-code the number of iterations
    const Int maxIter = 4; 

    Matrix<Real> scales, maxAbsValsB;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, scales );
        ColumnMaxNorms( B, maxAbsValsB );
        for( Int j=0; j<n; ++j )
            scales.Set
            ( j, 0, Max(scales.Get(j,0),maxAbsValsB.Get(j,0)) );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dCol );
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( RIGHT, NORMAL, scales, B );

        // Rescale the rows
        // ----------------
        RowMaxNorms( A, scales );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        RowMaxNorms( B, scales );
        ConeAllReduce( scales, orders, firstInds, mpi::MAX );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );
    }
    SetIndent( indent );
}

template<typename F>
void ConeRuizEquil
(       DistSparseMatrix<F>& A, 
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA, 
        DistMultiVec<Base<F>>& dRowB, 
        DistMultiVec<Base<F>>& dCol, 
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeRuizEquil"))
    typedef Base<F> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    dRowA.SetComm( comm );
    dRowB.SetComm( comm );
    dCol.SetComm( comm );
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose to control structure
    // For, simply hard-code a small number of iterations
    const Int maxIter = 4;

    DistMultiVec<Real> scales(comm), maxAbsValsB(comm);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, scales );
        ColumnMaxNorms( B, maxAbsValsB );
        for( Int j=0; j<n; ++j )
            scales.Set
            ( j, 0, Max(scales.Get(j,0),maxAbsValsB.Get(j,0)) );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dCol );
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( RIGHT, NORMAL, scales, B );

        // Rescale the rows 
        // ----------------
        RowMaxNorms( A, scales );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        RowMaxNorms( B, scales );
        ConeAllReduce( scales, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );
    }
    SetIndent( indent );
}

#define PROTO(F) \
  template void ConeRuizEquil \
  (       Matrix<F>& A, \
          Matrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void ConeRuizEquil \
  (       AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<Base<F>>& dRowA, \
          AbstractDistMatrix<Base<F>>& dRowB, \
          AbstractDistMatrix<Base<F>>& dCol, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff, bool progress ); \
  template void ConeRuizEquil \
  (       SparseMatrix<F>& A, \
          SparseMatrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void ConeRuizEquil \
  (       DistSparseMatrix<F>& A, \
          DistSparseMatrix<F>& B, \
          DistMultiVec<Base<F>>& dRowA, \
          DistMultiVec<Base<F>>& dRowB, \
          DistMultiVec<Base<F>>& dCol, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff, bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
