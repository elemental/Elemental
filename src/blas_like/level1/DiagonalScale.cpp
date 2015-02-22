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
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int ldim = A.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const T delta = d.Get(i,0);
            T* ABuffer = A.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    ABuffer[j*ldim] *= Conj(delta);
            else
                for( Int j=0; j<n; ++j )
                    ABuffer[j*ldim] *= delta;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const T delta = d.Get(j,0);
            T* ABuffer = A.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    ABuffer[i] *= Conj(delta);
            else
                for( Int i=0; i<m; ++i )
                    ABuffer[i] *= delta;
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<TDiag>& dPre, DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( side == LEFT )
    {
        ProxyCtrl ctrl;
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
        ProxyCtrl ctrl;
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
  const AbstractDistMatrix<TDiag>& d, AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
        DiagonalScale( side, orientation, d, ACast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const Matrix<TDiag>& d, SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const bool conjugate = ( orientation == ADJOINT );
    T* vBuf = A.ValueBuffer();
    const Int numEntries = A.NumEntries();
    if( side == LEFT )
    {
        if( d.Height() != A.Height() )
            LogicError("The size of d must match the height of A");
        for( Int k=0; k<numEntries; ++k )
        {
            const Int i = A.Row(k);
            const TDiag delta = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            vBuf[k] *= T(delta);
        }
    }
    else
    {
        if( d.Height() != A.Width() )
            LogicError("The size of d must match the width of A");
        for( Int k=0; k<numEntries; ++k )
        {
            const Int j = A.Col(k);
            const TDiag delta = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            vBuf[k] *= T(delta);
        }
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d, DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    if( !mpi::Congruent( d.Comm(), A.Comm() ) )
        LogicError("Communicators must be congruent");
    const bool conjugate = ( orientation == ADJOINT );
    if( side == LEFT )
    {
        if( d.Height() != A.Height() )
            LogicError("The size of d must match the height of A");
        // TODO: Ensure that the DistMultiVec conforms
        T* vBuf = A.ValueBuffer();
        const Int firstLocalRow = d.FirstLocalRow();
        for( Int k=0; k<A.NumLocalEntries(); ++k )
        {
            const Int i = A.Row(k);
            const Int iLoc = i - firstLocalRow;
            const TDiag delta = 
              ( conjugate ? Conj(d.GetLocal(iLoc,0)) : d.GetLocal(iLoc,0) );
            vBuf[k] *= T(delta);
        }
    }
    else
    {
        if( d.Height() != A.Width() )
            LogicError("The size of d must match the width of A");
        // NOTE: This is likely grossly suboptimal
        DistSparseMatrix<T> ATrans;
        Transpose( A, ATrans, conjugate );
        DiagonalScale( LEFT, NORMAL, d, ATrans );
        Transpose( ATrans, A, conjugate );
    }
}

template<typename TDiag,typename T>
void DiagonalScale
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<TDiag>& d, DistMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalScale"))
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
    for( Int iLoc=0; iLoc<d.LocalHeight(); ++iLoc ) 
    {
        const T delta = 
            ( conjugate ? Conj(d.GetLocal(iLoc,0)) : d.GetLocal(iLoc,0) );
        for( Int j=0; j<width; ++j )
            X.SetLocal( iLoc, j, delta*X.GetLocal(iLoc,j) ); 
    }
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<T,U,V>& A );

#define DIST_PROTO_REAL(T,U,V) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<Complex<T>,U,V>& A );

#define PROTO(T) \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, SparseMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistSparseMatrix<T>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistMultiVec<T>& X ); \
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
    const Matrix<T>& d, Matrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, SparseMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistSparseMatrix<Complex<T>>& A ); \
  template void DiagonalScale \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistMultiVec<Complex<T>>& X ); \
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

#include "El/macros/Instantiate.h"

} // namespace El
