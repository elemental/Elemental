/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, 
        Matrix<F>& A, 
  bool checkIfSingular )
{
    DEBUG_ONLY(CSE cse("DiagonalSolve"))
    const Int m = A.Height();
    const Int n = A.Width();
    const bool conj = ( orientation == ADJOINT );
    F* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    const FDiag* dBuf = d.LockedBuffer();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const F delta = ( conj ? Conj(dBuf[i]) : dBuf[i] );
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            for( Int j=0; j<n; ++j )
                ABuf[i+j*ALDim] *= deltaInv;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const F delta = ( conj ? Conj(dBuf[j]) : dBuf[j] );
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            for( Int i=0; i<m; ++i )
                ABuf[i+j*ALDim] *= deltaInv;
        }
    }
}

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& dPre, 
        DistMatrix<F,U,V>& A,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("DiagonalSolve");
      AssertSameGrids( dPre, A );
    )
    if( side == LEFT )
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.ColAlign();
        auto dPtr = ReadProxy<FDiag,U,Collect<V>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalSolve
        ( LEFT, orientation, d.LockedMatrix(), A.Matrix(), checkIfSingular );
    }
    else
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.RowAlign();
        auto dPtr = ReadProxy<FDiag,V,Collect<U>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalSolve
        ( RIGHT, orientation, d.LockedMatrix(), A.Matrix(), checkIfSingular );
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& d,
        AbstractDistMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_ONLY(CSE cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<DistMatrix<F,CDIST,RDIST>&>(A); \
        DiagonalSolve( side, orientation, d, ACast, checkIfSingular );
    #include "El/macros/GuardAndPayload.h"
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d,
        SparseMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("DiagonalSolve");
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

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d,
        DistSparseMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("DiagonalSolve");
      if( d.Width() != 1 )
          LogicError("d must be a column vector");
      if( !mpi::Congruent( d.Comm(), A.Comm() ) )
          LogicError("Communicators must be congruent");
    )
    const bool conjugate = ( orientation == ADJOINT );
    F* vBuf = A.ValueBuffer();
    const Int numEntries = A.NumLocalEntries();
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
            const Int i = A.Row(k);
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
        vector<F> sendVals( numSendInds );
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            sendVals[s] = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
        }

        // Now send them
        vector<F> recvVals( meta.numRecvInds );
        mpi::AllToAll
        ( sendVals.data(), meta.sendSizes.data(), meta.sendOffs.data(),
          recvVals.data(), meta.recvSizes.data(), meta.recvOffs.data(), 
          A.Comm() );

        // Loop over the entries of A and rescale
        for( Int k=0; k<numEntries; ++k )
            vBuf[k] /= recvVals[meta.colOffs[k]];
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d,
        DistMultiVec<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CSE cse("DiagonalSolve"))
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
    const FDiag* dBuf = d.LockedMatrix().LockedBuffer();
    F* XBuf = X.Matrix().Buffer();
    const Int XLDim = X.Matrix().LDim();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const F delta = ( conjugate ? Conj(dBuf[iLoc]) : dBuf[iLoc] );
        DEBUG_ONLY(
          if( checkIfSingular && delta == F(0) )
              throw SingularMatrixException(); 
        )
        for( Int j=0; j<width; ++j )
            XBuf[iLoc+j*XLDim] /= delta;
    }
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, \
          DistMatrix<T,U,V>& A, \
    bool checkIfSingular );

#define DIST_PROTO_REAL(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, \
          DistMatrix<Complex<T>,U,V>& A, \
    bool checkIfSingular );

#define PROTO(T) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<T>& A, \
     bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, \
          AbstractDistMatrix<T>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<T>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<T>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistMultiVec<T>& X, \
    bool checkIfSingular ); \
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
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, \
          AbstractDistMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistMultiVec<Complex<T>>& X, \
    bool checkIfSingular ); \
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

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
