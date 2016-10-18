/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real DampScaling( Real alpha )
{
    const Real tol = Pow(limits::Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else 
        return Max(alpha,tol);
}

namespace cone {

template<typename F>
void RuizEquil
(       Matrix<F>& A, 
        Matrix<F>& B, 
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB, 
        Matrix<Base<F>>& dCol, 
  const Matrix<Int>& orders,  
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_CSE
    LogicError("This routine is not yet written");
}

template<typename F>
void RuizEquil
(       ElementalMatrix<F>& APre, 
        ElementalMatrix<F>& BPre,
        ElementalMatrix<Base<F>>& dRowAPre, 
        ElementalMatrix<Base<F>>& dRowBPre,
        ElementalMatrix<Base<F>>& dColPre,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    DEBUG_CSE
    typedef Base<F> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<F,F,MC,MR>
      AProx( APre, control ),
      BProx( BPre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR>
      dRowAProx( dRowAPre, control ),
      dRowBProx( dRowBPre, control );
    DistMatrixWriteProxy<Real,Real,MR,STAR>
      dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& dRowA = dRowAProx.Get();
    auto& dRowB = dRowBProx.Get();
    auto& dCol = dColProx.Get();

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
    auto& colScaleLoc = colScale.Matrix();
    auto& colScaleBLoc = colScaleB.Matrix();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Rescale the columns
        // -------------------
        ColumnMaxNorms( A, colScale );
        ColumnMaxNorms( B, colScaleB );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            colScaleLoc(jLoc) = Max(colScaleLoc(jLoc),colScaleBLoc(jLoc));
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
        cone::AllReduce( rowScale, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( rowScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScale, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScale, B );
    }
    SetIndent( indent );
}

template<typename F>
void RuizEquil
(       SparseMatrix<F>& A, 
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_CSE
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
            scales(j) = Max(scales(j),maxAbsValsB(j));
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
        cone::AllReduce( scales, orders, firstInds, mpi::MAX );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );
    }
    SetIndent( indent );
}

template<typename F>
void RuizEquil
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
    DEBUG_CSE
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

        Real* scaleBuf = scales.Matrix().Buffer();
        const Real* maxAbsValBuf = maxAbsValsB.LockedMatrix().LockedBuffer();
        const Int nLoc = scales.LocalHeight();
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            scaleBuf[jLoc] = Max(scaleBuf[jLoc],maxAbsValBuf[jLoc]);

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
        cone::AllReduce( scales, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );
    }
    SetIndent( indent );
}

#define PROTO(F) \
  template void RuizEquil \
  (       Matrix<F>& A, \
          Matrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void RuizEquil \
  (       ElementalMatrix<F>& A, \
          ElementalMatrix<F>& B, \
          ElementalMatrix<Base<F>>& dRowA, \
          ElementalMatrix<Base<F>>& dRowB, \
          ElementalMatrix<Base<F>>& dCol, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
    Int cutoff, bool progress ); \
  template void RuizEquil \
  (       SparseMatrix<F>& A, \
          SparseMatrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void RuizEquil \
  (       DistSparseMatrix<F>& A, \
          DistSparseMatrix<F>& B, \
          DistMultiVec<Base<F>>& dRowA, \
          DistMultiVec<Base<F>>& dRowB, \
          DistMultiVec<Base<F>>& dCol, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff, bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace cone
} // namespace El
