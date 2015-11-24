/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
    const Int m = A.Height();
    const Int n = A.Width();
    const bool conj = ( orientation == ADJOINT );
    const TDiag* dBuf = d.LockedBuffer();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const T delta = ( conj ? Conj(dBuf[i]) : dBuf[i] );
            for( Int j=0; j<n; ++j )
                ABuf[i+j*ALDim] *= delta;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T delta = ( conj ? Conj(dBuf[j]) : dBuf[j] );
            for( Int i=0; i<m; ++i )
                ABuf[i+j*ALDim] *= delta;
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<TDiag>& dPre, DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
    if( side == LEFT )
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.ColAlign();
        auto dPtr = ReadProxy<TDiag,U,Collect<V>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalScale( LEFT, orientation, d.LockedMatrix(), A.Matrix() );
    }
    else
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.RowAlign();
        auto dPtr = ReadProxy<TDiag,V,Collect<U>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalScale( RIGHT, orientation, d.LockedMatrix(), A.Matrix() );
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<TDiag>& d, ElementalMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
        DiagonalScale( side, orientation, d, ACast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, 
        SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
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
        if( d.Height() != A.Height() )
            LogicError("The size of d must match the height of A");
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = rowBuf[k];
            const T delta = ( conjugate ? Conj(dBuf[i]) : dBuf[i] );
            vBuf[k] *= delta;
        }
    }
    else
    {
        if( d.Height() != A.Width() )
            LogicError("The size of d must match the width of A");
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
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d,
        DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    if( !mpi::Congruent( d.Comm(), A.Comm() ) )
        LogicError("Communicators must be congruent");
    const bool conjugate = ( orientation == ADJOINT );
    const Int numEntries = A.NumLocalEntries();
    T* vBuf = A.ValueBuffer();
    const Int* rowBuf = A.LockedSourceBuffer();
    const TDiag* dBuf = d.LockedMatrix().LockedBuffer();
    const Int firstLocalRow = d.FirstLocalRow();
    if( side == LEFT )
    {
        DEBUG_ONLY(
          if( d.Height() != A.Height() )
              LogicError("The size of d must match the height of A");
        )
        // TODO: Ensure that the DistMultiVec conforms
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
        DEBUG_ONLY(
          if( d.Height() != A.Width() )
              LogicError("The size of d must match the width of A");
        )
        A.InitializeMultMeta();
        const auto& meta = A.multMeta;

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
          A.Comm() );

        // Loop over the entries of A and rescale
        for( Int k=0; k<numEntries; ++k )
            vBuf[k] *= recvVals[meta.colOffs[k]];
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d,
        DistMultiVec<T>& X )
{
    DEBUG_ONLY(
      CSE cse("DiagonalScale");
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( !mpi::Congruent( d.Comm(), X.Comm() ) )
          LogicError("Communicators must be congruent");
      if( side != LEFT )
          LogicError("Only the 'LEFT' argument is currently supported");
      if( d.Height() != X.Height() )
          LogicError("d and X must be the same size");
    )
    const bool conjugate = ( orientation == ADJOINT );
    const Int width = X.Width();
    T* XBuf = X.Matrix().Buffer();
    const Int XLDim = X.Matrix().LDim();
    const Int localHeight = d.LocalHeight();
    const TDiag* dBuf = d.LockedMatrix().LockedBuffer();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
    {
        const T delta = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
        for( Int j=0; j<width; ++j )
            XBuf[iLoc+j*XLDim] *= delta;
    }
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<T>& d, \
          DistMatrix<T,U,V>& A );

#define DIST_PROTO_REAL(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<T>& d, \
          DistMatrix<Complex<T>,U,V>& A );

#define PROTO(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<T>& d, \
          ElementalMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistMultiVec<T>& X ); \
  DIST_PROTO(T,CIRC,CIRC); \
  DIST_PROTO(T,MC,  MR  ); \
  DIST_PROTO(T,MC,  STAR); \
  DIST_PROTO(T,MD,  STAR); \
  DIST_PROTO(T,MR,  MC  ); \
  DIST_PROTO(T,MR,  STAR); \
  DIST_PROTO(T,STAR,MC  ); \
  DIST_PROTO(T,STAR,MD  ); \
  DIST_PROTO(T,STAR,MR  ); \
  DIST_PROTO(T,STAR,STAR); \
  DIST_PROTO(T,STAR,VC  ); \
  DIST_PROTO(T,STAR,VR  ); \
  DIST_PROTO(T,VC  ,STAR); \
  DIST_PROTO(T,VR  ,STAR);

#define PROTO_REAL(T) \
  PROTO(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const ElementalMatrix<T>& d, \
          ElementalMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistMultiVec<Complex<T>>& X ); \
  DIST_PROTO_REAL(T,CIRC,CIRC); \
  DIST_PROTO_REAL(T,MC,  MR  ); \
  DIST_PROTO_REAL(T,MC,  STAR); \
  DIST_PROTO_REAL(T,MD,  STAR); \
  DIST_PROTO_REAL(T,MR,  MC  ); \
  DIST_PROTO_REAL(T,MR,  STAR); \
  DIST_PROTO_REAL(T,STAR,MC  ); \
  DIST_PROTO_REAL(T,STAR,MD  ); \
  DIST_PROTO_REAL(T,STAR,MR  ); \
  DIST_PROTO_REAL(T,STAR,STAR); \
  DIST_PROTO_REAL(T,STAR,VC  ); \
  DIST_PROTO_REAL(T,STAR,VR  ); \
  DIST_PROTO_REAL(T,VC  ,STAR); \
  DIST_PROTO_REAL(T,VR  ,STAR);

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
