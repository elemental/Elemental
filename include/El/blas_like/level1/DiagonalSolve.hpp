/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_DIAGONALSOLVE_HPP
#define EL_BLAS_DIAGONALSOLVE_HPP

namespace El {

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side,
  Orientation orientation,
  const Matrix<FDiag>& d, 
        Matrix<F>& A, 
  bool checkIfSingular )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const bool conj = ( orientation == ADJOINT );
    if( side == LEFT )
    {
        DEBUG_ONLY(
          if( d.Height() != m )
              LogicError("Invalid left diagonal solve dimension");
        )
        for( Int i=0; i<m; ++i )
        {
            const F delta = ( conj ? Conj(d(i)) : d(i) );
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            for( Int j=0; j<n; ++j )
                A(i,j) *= deltaInv;
        }
    }
    else
    {
        DEBUG_ONLY(
          if( d.Height() != n )
              LogicError("Invalid right diagonal solve dimension");
        )
        for( Int j=0; j<n; ++j )
        {
            const F delta = ( conj ? Conj(d(j)) : d(j) );
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            for( Int i=0; i<m; ++i )
                A(i,j) *= deltaInv;
        }
    }
}

template<typename F>
void SymmetricDiagonalSolve
( const Matrix<Base<F>>& d, 
        Matrix<F>& A )
{
    DEBUG_CSE
    const Int n = A.Width();
    DEBUG_ONLY(
      if( d.Height() != n )
          LogicError("Invalid symmetric diagonal solve dimension");
    )
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i ) 
            A(i,j) /= d(i)*d(j);
}

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<FDiag>& dPre, 
        DistMatrix<F,U,V>& A,
  bool checkIfSingular )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( dPre, A );
    )
    if( side == LEFT )
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.ColAlign();

        DistMatrixReadProxy<FDiag,FDiag,U,Collect<V>()> dProx( dPre, ctrl );
        auto& d = dProx.GetLocked();

        DiagonalSolve
        ( LEFT, orientation, d.LockedMatrix(), A.Matrix(), checkIfSingular );
    }
    else
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.RowAlign();

        DistMatrixReadProxy<FDiag,FDiag,V,Collect<U>()> dProx( dPre, ctrl );
        auto& d = dProx.GetLocked();

        DiagonalSolve
        ( RIGHT, orientation, d.LockedMatrix(), A.Matrix(), checkIfSingular );
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<FDiag>& d,
        ElementalMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_CSE
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = static_cast<DistMatrix<F,CDIST,RDIST>&>(A); \
        DiagonalSolve( side, orientation, d, ACast, checkIfSingular );
    #include <El/macros/GuardAndPayload.h>
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d,
        SparseMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
    )
    const bool conjugate = ( orientation == ADJOINT );
    const Int numEntries = A.NumEntries();
    F* vBuf = A.ValueBuffer();
    const Int* rowBuf = A.LockedSourceBuffer();
    const Int* colBuf = A.LockedTargetBuffer();
    const FDiag* dBuf = d.LockedBuffer();
    if( side == LEFT )
    {
        DEBUG_ONLY(
          if( d.Height() != A.Height() )
              LogicError("The size of d must match the height of A");
        )
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = rowBuf[k];
            const FDiag delta = ( conjugate ? Conj(dBuf[i]) : dBuf[i] );
            if( checkIfSingular && delta == FDiag(0) )
                throw SingularMatrixException();
            vBuf[k] /= F(delta);
        }
    }
    else
    {
        DEBUG_ONLY(
          if( d.Height() != A.Width() )
              LogicError("The size of d must match the width of A");
        )
        for( Int k=0; k<numEntries; ++k )
        {
            const Int j = colBuf[k];
            const FDiag delta = ( conjugate ? Conj(dBuf[j]) : dBuf[j] );
            if( checkIfSingular && delta == FDiag(0) )
                throw SingularMatrixException();
            vBuf[k] /= F(delta);
        }
    }
}

template<typename F>
void SymmetricDiagonalSolve( const Matrix<Base<F>>& d, SparseMatrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( d.Height() != A.Height() )
          LogicError("The size of d must match the height of A");
    )
    typedef Base<F> Real;

    const Int numEntries = A.NumEntries();
    F* vBuf = A.ValueBuffer();
    const Int* rowBuf = A.LockedSourceBuffer();
    const Int* colBuf = A.LockedTargetBuffer();

    const Real* dBuf = d.LockedBuffer();

    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = rowBuf[k];
        const Int j = colBuf[k];
        const Real deltaRow = dBuf[i];
        const Real deltaCol = dBuf[j];
        DEBUG_ONLY(
          if( deltaRow*deltaCol == Real(0) )
              throw SingularMatrixException();
        )
        vBuf[k] /= deltaRow*deltaCol;
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d,
        DistSparseMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( !mpi::Congruent( d.Comm(), A.Comm() ) )
          LogicError("Communicators must be congruent");
    )
    const bool conjugate = ( orientation == ADJOINT );

    const Int numEntries = A.NumLocalEntries();
    F* vBuf = A.ValueBuffer();
    const Int* rBuf = A.LockedSourceBuffer();

    const FDiag* dBuf = d.LockedMatrix().LockedBuffer();
    const Int firstLocalRow = d.FirstLocalRow();

    if( side == LEFT )
    {
        DEBUG_ONLY(
          if( d.Height() != A.Height() )
              LogicError("The length of d must match the height of A");
        )
        // TODO: Ensure that the DistMultiVec conforms
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = rBuf[k];
            const Int iLoc = i - firstLocalRow;
            const F delta = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
            DEBUG_ONLY(
              if( checkIfSingular && delta == F(0) )
                  throw SingularMatrixException();
            )
            vBuf[k] /= delta;
        }
    }
    else
    {
        A.InitializeMultMeta();
        const auto& meta = A.multMeta;

        // Pack the send values
        const Int numSendInds = meta.sendInds.size();
        vector<F> sendVals;
        FastResize( sendVals, numSendInds );
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            sendVals[s] = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
        }

        // Now send them
        vector<F> recvVals;
        FastResize( recvVals, meta.numRecvInds );
        mpi::AllToAll
        ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
          recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(), 
          A.Comm() );

        // Loop over the entries of A and rescale
        for( Int k=0; k<numEntries; ++k )
            vBuf[k] /= recvVals[meta.colOffs[k]];
    }
}

template<typename F>
void SymmetricDiagonalSolve
( const DistMultiVec<Base<F>>& d,
        DistSparseMatrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( d.Height() != A.Height() )
          LogicError("The length of d must match the height of A");
      if( !mpi::Congruent( d.Comm(), A.Comm() ) )
          LogicError("Communicators must be congruent");
    )
    typedef Base<F> Real;

    const Int numEntries = A.NumLocalEntries();
    F* vBuf = A.ValueBuffer();
    const Int* rBuf = A.LockedSourceBuffer();

    const Real* dBuf = d.LockedMatrix().LockedBuffer();
    const Int firstLocalRow = d.FirstLocalRow();

    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

    // Pack the send values
    const Int numSendInds = meta.sendInds.size();
    vector<Real> sendVals( numSendInds );
    for( Int s=0; s<numSendInds; ++s )
    {
        const Int i = meta.sendInds[s];
        const Int iLoc = i - firstLocalRow;
        sendVals[s] = dBuf[iLoc];
    }

    // Now send them
    vector<Real> recvVals( meta.numRecvInds );
    mpi::AllToAll
    ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
      recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(), 
      A.Comm() );

    // Loop over the entries of A and rescale
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = rBuf[k];
        const Int iLoc = i - firstLocalRow;
        vBuf[k] /= recvVals[meta.colOffs[k]]*dBuf[iLoc];
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d,
        DistMultiVec<F>& X,
  bool checkIfSingular )
{
    DEBUG_CSE
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    if( !mpi::Congruent( d.Comm(), X.Comm() ) )
        LogicError("Communicators must be congruent");
    if( side != LEFT )
        LogicError("Only the 'LEFT' argument is currently supported");
    if( d.Height() != X.Height() )
        LogicError("d and X must be the same size");
    const bool conjugate = ( orientation == ADJOINT );
    const Int width = X.Width();
    const Int localHeight = d.LocalHeight();
    auto& XLoc = X.Matrix();
    auto& dLoc = d.LockedMatrix();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const F delta = ( conjugate ? Conj(dLoc(iLoc)) : dLoc(iLoc) );
        DEBUG_ONLY(
          if( checkIfSingular && delta == F(0) )
              throw SingularMatrixException(); 
        )
        for( Int j=0; j<width; ++j )
            XLoc(iLoc,j) /= delta;
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(F) \
  EL_EXTERN template void DiagonalSolve \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<F>& d, \
          Matrix<F>& A, \
    bool checkIfSingular ); \
  EL_EXTERN template void SymmetricDiagonalSolve \
  ( const Matrix<Base<F>>& d, \
          Matrix<F>& A ); \
  EL_EXTERN template void DiagonalSolve \
  ( LeftOrRight side, \
    Orientation orientation, \
    const ElementalMatrix<F>& d, \
          ElementalMatrix<F>& A, \
    bool checkIfSingular ); \
  EL_EXTERN template void DiagonalSolve \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<F>& d, \
          SparseMatrix<F>& A, \
    bool checkIfSingular ); \
  EL_EXTERN template void SymmetricDiagonalSolve \
  ( const Matrix<Base<F>>& d, SparseMatrix<F>& A ); \
  EL_EXTERN template void DiagonalSolve \
  ( LeftOrRight side, \
    Orientation orientation, \
    const DistMultiVec<F>& d, \
          DistSparseMatrix<F>& A, \
    bool checkIfSingular ); \
  EL_EXTERN template void SymmetricDiagonalSolve \
  ( const DistMultiVec<Base<F>>& d, \
          DistSparseMatrix<F>& A ); \
  EL_EXTERN template void DiagonalSolve \
  ( LeftOrRight side, \
    Orientation orientation, \
    const DistMultiVec<F>& d, \
          DistMultiVec<F>& X, \
    bool checkIfSingular );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_DIAGONALSOLVE_HPP
