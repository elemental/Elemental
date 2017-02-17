/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_DIAGONALSCALE_HPP
#define EL_BLAS_DIAGONALSCALE_HPP

namespace El {

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const Matrix<TDiag>& d,
        Matrix<T>& A )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const bool conj = ( orientation == ADJOINT );
    if( side == LEFT )
    {
        EL_DEBUG_ONLY(
          if( d.Height() != m )
              LogicError("Invalid left diagonal scaling dimension");
        )
        for( Int i=0; i<m; ++i )
        {
            const T delta = ( conj ? Conj(d(i)) : d(i) );
            for( Int j=0; j<n; ++j )
                A(i,j) *= delta;
        }
    }
    else
    {
        EL_DEBUG_ONLY(
          if( d.Height() != n )
              LogicError("Invalid right diagonal scaling dimension");
        )
        for( Int j=0; j<n; ++j )
        {
            const T delta = ( conj ? Conj(d(j)) : d(j) );
            for( Int i=0; i<m; ++i )
                A(i,j) *= delta;
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V,DistWrap wrapType>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const AbstractDistMatrix<TDiag>& dPre,
        DistMatrix<T,U,V,wrapType>& A )
{
    EL_DEBUG_CSE
    if( wrapType == ELEMENT )
    {
        if( side == LEFT )
        {
            ElementalProxyCtrl ctrl;
            ctrl.rootConstrain = true;
            ctrl.colConstrain = true;
            ctrl.root = A.Root();
            ctrl.colAlign = A.ColAlign();

            DistMatrixReadProxy<TDiag,TDiag,U,Collect<V>()> dProx( dPre, ctrl );
            auto& d = dProx.GetLocked();

            DiagonalScale( LEFT, orientation, d.LockedMatrix(), A.Matrix() );
        }
        else
        {
            ElementalProxyCtrl ctrl;
            ctrl.rootConstrain = true;
            ctrl.colConstrain = true;
            ctrl.root = A.Root();
            ctrl.colAlign = A.RowAlign();

            DistMatrixReadProxy<TDiag,TDiag,V,Collect<U>()> dProx( dPre, ctrl );
            auto& d = dProx.GetLocked();

            DiagonalScale( RIGHT, orientation, d.LockedMatrix(), A.Matrix() );
        }
    }
    else
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();

        if( side == LEFT )
        {
            ctrl.colAlign = A.ColAlign();
            ctrl.blockHeight = A.BlockHeight();
            ctrl.colCut = A.ColCut();

            DistMatrixReadProxy<TDiag,TDiag,U,Collect<V>(),BLOCK>
              dProx( dPre, ctrl );
            auto& d = dProx.GetLocked();

            DiagonalScale( LEFT, orientation, d.LockedMatrix(), A.Matrix() );
        }
        else
        {
            ctrl.colAlign = A.RowAlign();
            ctrl.blockHeight = A.BlockWidth();
            ctrl.colCut = A.RowCut();

            DistMatrixReadProxy<TDiag,TDiag,V,Collect<U>(),BLOCK>
              dProx( dPre, ctrl );
            auto& d = dProx.GetLocked();

            DiagonalScale( RIGHT, orientation, d.LockedMatrix(), A.Matrix() );
        }
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const AbstractDistMatrix<TDiag>& d,
        AbstractDistMatrix<T>& A )
{
    EL_DEBUG_CSE
    #define GUARD(CDIST,RDIST,WRAP) \
      A.ColDist() == CDIST && A.RowDist() == RDIST && A.Wrap() == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
        auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST,WRAP>&>(A); \
        DiagonalScale( side, orientation, d, ACast );
    #include <El/macros/GuardAndPayload.h>
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const Matrix<TDiag>& d,
        SparseMatrix<T>& A )
{
    EL_DEBUG_CSE
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const bool conjugate = ( orientation == ADJOINT );
    const Int numEntries = A.NumEntries();
    T* vBuf = A.ValueBuffer();
    const Int* rowBuf = A.LockedSourceBuffer();
    const Int* colBuf = A.LockedTargetBuffer();
    const TDiag* dBuf = d.LockedBuffer();
    if( side == LEFT )
    {
        EL_DEBUG_ONLY(
          if( d.Height() != A.Height() )
              LogicError("The size of d must match the height of A");
        )
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = rowBuf[k];
            const T delta = ( conjugate ? Conj(dBuf[i]) : dBuf[i] );
            vBuf[k] *= delta;
        }
    }
    else
    {
        EL_DEBUG_ONLY(
          if( d.Height() != A.Width() )
              LogicError("The size of d must match the width of A");
        )
        for( Int k=0; k<numEntries; ++k )
        {
            const Int j = colBuf[k];
            const T delta = ( conjugate ? Conj(dBuf[j]) : dBuf[j] );
            vBuf[k] *= delta;
        }
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const DistMultiVec<TDiag>& d,
        DistSparseMatrix<T>& A )
{
    EL_DEBUG_CSE
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    if( !mpi::Congruent( d.Grid().Comm(), A.Grid().Comm() ) )
        LogicError("Communicators must be congruent");
    const bool conjugate = ( orientation == ADJOINT );
    const Int numEntries = A.NumLocalEntries();
    T* vBuf = A.ValueBuffer();
    const Int* rowBuf = A.LockedSourceBuffer();
    const TDiag* dBuf = d.LockedMatrix().LockedBuffer();
    const Int firstLocalRow = d.FirstLocalRow();
    if( side == LEFT )
    {
        EL_DEBUG_ONLY(
          if( d.Height() != A.Height() )
              LogicError("The size of d must match the height of A");
        )
        // TODO(poulson): Ensure that the DistMultiVec conforms
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = rowBuf[k];
            const Int iLoc = i - firstLocalRow;
            const T delta = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
            vBuf[k] *= delta;
        }
    }
    else
    {
        EL_DEBUG_ONLY(
          if( d.Height() != A.Width() )
              LogicError("The size of d must match the width of A");
        )
        A.InitializeMultMeta();
        const auto& meta = A.LockedDistGraph().multMeta;
        // Pack the send values
        const Int numSendInds = meta.sendInds.size();
        vector<T> sendVals;
        FastResize( sendVals, numSendInds );

        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            sendVals[s] = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
        }

        // Now send them
        vector<T> recvVals;
        FastResize( recvVals, meta.numRecvInds );
        mpi::AllToAll
        ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
          recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(),
          A.Grid().Comm() );

        // Loop over the entries of A and rescale
        for( Int k=0; k<numEntries; ++k )
            vBuf[k] *= recvVals[meta.colOffs[k]];
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side,
  Orientation orientation,
  const DistMultiVec<TDiag>& d,
        DistMultiVec<T>& X )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( !mpi::Congruent( d.Grid().Comm(), X.Grid().Comm() ) )
          LogicError("Communicators must be congruent");
      if( side != LEFT )
          LogicError("Only the 'LEFT' argument is currently supported");
      if( d.Height() != X.Height() )
          LogicError("d and X must be the same size");
    )
    const bool conjugate = ( orientation == ADJOINT );
    const Int width = X.Width();
    auto& XLoc = X.Matrix();
    auto& dLoc = d.LockedMatrix();
    const Int localHeight = d.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const T delta = ( conjugate ? Conj(dLoc(iLoc)) : dLoc(iLoc) );
        for( Int j=0; j<width; ++j )
            XLoc(iLoc,j) *= delta;
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void DiagonalScale \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<T>& A ); \
  EL_EXTERN template void DiagonalScale \
  ( LeftOrRight side, \
    Orientation orientation, \
    const AbstractDistMatrix<T>& d, \
          AbstractDistMatrix<T>& A ); \
  EL_EXTERN template void DiagonalScale \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<T>& A ); \
  EL_EXTERN template void DiagonalScale \
  ( LeftOrRight side, \
    Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<T>& A ); \
  EL_EXTERN template void DiagonalScale \
  ( LeftOrRight side, \
    Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistMultiVec<T>& X );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_DIAGONALSCALE_HPP
