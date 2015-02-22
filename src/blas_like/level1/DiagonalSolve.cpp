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
  const Matrix<FDiag>& d, Matrix<F>& A, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int ldim = A.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* ABuffer = A.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    ABuffer[j*ldim] *= Conj(deltaInv);
            else
                for( Int j=0; j<n; ++j )
                    ABuffer[j*ldim] *= deltaInv;
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const F delta = d.Get(j,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* ABuffer = A.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    ABuffer[i] *= Conj(deltaInv);
            else
                for( Int i=0; i<m; ++i )
                    ABuffer[i] *= deltaInv;
        }
    }
}

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& dPre, 
        DistMatrix<F,U,V>& A, bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("DiagonalSolve");
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
  const AbstractDistMatrix<FDiag>& d, AbstractDistMatrix<F>& A,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = dynamic_cast<DistMatrix<F,CDIST,RDIST>&>(A); \
        DiagonalSolve( side, orientation, d, ACast, checkIfSingular );
    #include "El/macros/GuardAndPayload.h"
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const Matrix<FDiag>& d, SparseMatrix<F>& A, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    const bool conjugate = ( orientation == ADJOINT );
    F* vBuf = A.ValueBuffer();
    if( side == LEFT )
    {
        if( d.Height() != A.Height() )
            LogicError("The size of d must match the height of A");
        for( Int k=0; k<A.NumEntries(); ++k )
        {
            const Int i = A.Row(k);
            const FDiag delta = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            if( checkIfSingular && delta == FDiag(0) )
                throw SingularMatrixException();
            vBuf[k] /= F(delta);
        }
    }
    else
    {
        if( d.Height() != A.Width() )
            LogicError("The size of d must match the width of A");
        for( Int k=0; k<A.NumEntries(); ++k )
        {
            const Int j = A.Col(k);
            const FDiag delta = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            if( checkIfSingular && delta == FDiag(0) )
                throw SingularMatrixException();
            vBuf[k] /= F(delta);
        }
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d, DistSparseMatrix<F>& A, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    if( d.Width() != 1 )
        LogicError("d must be a column vector");
    if( !mpi::Congruent( d.Comm(), A.Comm() ) )
        LogicError("Communicators must be congruent");
    const bool conjugate = ( orientation == ADJOINT );
    if( side == LEFT )
    {
        if( d.Height() != A.Height() )
            LogicError("The length of d must match the height of A");
        // TODO: Ensure that the DistMultiVec conforms
        F* vBuf = A.ValueBuffer();
        const Int firstLocalRow = d.FirstLocalRow();
        for( Int k=0; k<A.NumLocalEntries(); ++k )
        {
            const Int i = A.Row(k);
            const Int iLoc = i - firstLocalRow;
            const FDiag delta = 
              ( conjugate ? Conj(d.GetLocal(iLoc,0)) : d.GetLocal(iLoc,0) );
            if( checkIfSingular && delta == FDiag(0) )
                throw SingularMatrixException();
            vBuf[k] /= F(delta);
        }
    }
    else
    {
        // NOTE: This is likely grossly suboptimal
        DistSparseMatrix<F> ATrans;
        Transpose( A, ATrans, conjugate );
        DiagonalSolve( LEFT, NORMAL, d, ATrans );
        Transpose( ATrans, A, conjugate );
    }
}

template<typename FDiag,typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMultiVec<FDiag>& d, DistMultiVec<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
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
        const F delta = 
            ( conjugate ? Conj(d.GetLocal(iLoc,0)) : d.GetLocal(iLoc,0) );
        if( checkIfSingular && delta == F(0) )
            throw SingularMatrixException(); 
        for( Int j=0; j<width; ++j )
            X.SetLocal( iLoc, j, X.GetLocal(iLoc,j)/delta );
    }
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<T,U,V>& A, \
    bool checkIfSingular );

#define DIST_PROTO_REAL(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<Complex<T>,U,V>& A, \
    bool checkIfSingular );

#define PROTO(T) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& A, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, SparseMatrix<T>& A, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistSparseMatrix<T>& A, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistMultiVec<T>& X, bool checkIfSingular ); \
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
    const Matrix<T>& d, Matrix<Complex<T>>& A, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, SparseMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistSparseMatrix<Complex<T>>& A, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMultiVec<T>& d, DistMultiVec<Complex<T>>& X, \
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
#include "El/macros/Instantiate.h"

} // namespace El
