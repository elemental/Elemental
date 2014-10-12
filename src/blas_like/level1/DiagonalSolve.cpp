/*
   Copyright (c) 2009-2014, Jack Poulson
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
  const Matrix<FDiag>& d, Matrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int ldim = X.LDim();
    if( side == LEFT )
    {
        for( Int i=0; i<m; ++i )
        {
            const F delta = d.Get(i,0);
            if( checkIfSingular && delta == F(0) )
                throw SingularMatrixException();
            const F deltaInv = F(1)/delta;
            F* XBuffer = X.Buffer(i,0);
            if( orientation == ADJOINT )
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= Conj(deltaInv);
            else
                for( Int j=0; j<n; ++j )
                    XBuffer[j*ldim] *= deltaInv;
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
            F* XBuffer = X.Buffer(0,j);
            if( orientation == ADJOINT )
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= Conj(deltaInv);
            else
                for( Int i=0; i<m; ++i )
                    XBuffer[i] *= deltaInv;
        }
    }
}

template<typename FDiag,typename F,Dist U,Dist V>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<FDiag>& dPre, 
        DistMatrix<F,U,V>& X, bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("DiagonalSolve");
        AssertSameGrids( dPre, X );
    )
    if( side == LEFT )
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = X.Root();
        ctrl.colAlign = X.ColAlign();
        auto dPtr = ReadProxy<FDiag,U,GatheredDist<V>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalSolve
        ( LEFT, orientation, d.LockedMatrix(), X.Matrix(), checkIfSingular );
    }
    else
    {
        ProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = X.Root();
        ctrl.colAlign = X.RowAlign();
        auto dPtr = ReadProxy<FDiag,V,GatheredDist<U>()>( &dPre, ctrl );
        auto& d = *dPtr;
        DiagonalSolve
        ( RIGHT, orientation, d.LockedMatrix(), X.Matrix(), checkIfSingular );
    }
}

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) X.ColDist() == CDIST && X.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<F,CDIST,RDIST>&>(X); \
        DiagonalSolve( side, orientation, d, XCast, checkIfSingular );
    #include "El/macros/GuardAndPayload.h"
}

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<Complex<F>>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) X.ColDist() == CDIST && X.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<Complex<F>,CDIST,RDIST>&>(X); \
        DiagonalSolve( side, orientation, d, XCast, checkIfSingular );
    #include "El/macros/GuardAndPayload.h"
}

#define DIST_PROTO(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<T,U,V>& X, \
    bool checkIfSingular );

#define DIST_PROTO_REAL(T,U,V) \
  DIST_PROTO(T,U,V) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, DistMatrix<Complex<T>,U,V>& X, \
    bool checkIfSingular );

#define PROTO(T) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, \
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
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<Complex<T>>& X, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, \
    bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<Complex<T>>& X, \
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
