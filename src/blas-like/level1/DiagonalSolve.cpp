/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename FDiag,typename F>
inline void
DiagonalSolve
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

template<typename FDiag,typename F,Dist U,Dist V,Dist W,Dist Z>
inline void
DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const DistMatrix<FDiag,U,V>& d, DistMatrix<F,W,Z>& X, bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    if( side == LEFT )
    {
        const Dist ZGath = GatheredDist<Z>();
        if( U == W && V == STAR && d.ColAlign() == X.ColAlign() )
        {
            DiagonalSolve
            ( LEFT, orientation, d.LockedMatrix(), X.Matrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<FDiag,W,ZGath> d_W_ZGath( X.Grid() );
            d_W_ZGath = d;
            DiagonalSolve
            ( LEFT, orientation,
              d_W_ZGath.LockedMatrix(), X.Matrix(), checkIfSingular );
        }
    }
    else
    {
        const Dist WGath = GatheredDist<W>();
        if( U == Z && V == STAR && d.ColAlign() == X.RowAlign() )
        {
            DiagonalSolve
            ( RIGHT, orientation, d.LockedMatrix(), X.Matrix(),
              checkIfSingular );
        }
        else
        {
            DistMatrix<FDiag,Z,WGath> d_Z_WGath( X.Grid() );
            d_Z_WGath = d;
            DiagonalSolve
            ( RIGHT, orientation,
              d_Z_WGath.LockedMatrix(), X.Matrix(), checkIfSingular );
        }
    }
}

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) \
        d.DistData().colDist == CDIST && d.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& dCast = dynamic_cast<const DistMatrix<F,CDIST,RDIST>&>(d);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<F,CDIST,RDIST>&>(X); \
        DiagonalSolve( side, orientation, dCast, XCast, checkIfSingular );
    #include "El/core/NestedGuardAndPayload.h"
}

template<typename F>
void DiagonalSolve
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<F>& d, AbstractDistMatrix<Complex<F>>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(CallStackEntry cse("DiagonalSolve"))
    #define GUARD(CDIST,RDIST) \
        d.DistData().colDist == CDIST && d.DistData().rowDist == RDIST
    #define INNER_GUARD(CDIST,RDIST) \
        X.DistData().colDist == CDIST && X.DistData().rowDist == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& dCast = dynamic_cast<const DistMatrix<F,CDIST,RDIST>&>(d);
    #define INNER_PAYLOAD(CDIST,RDIST) \
        auto& XCast = dynamic_cast<DistMatrix<Complex<F>,CDIST,RDIST>&>(X); \
        DiagonalSolve( side, orientation, dCast, XCast, checkIfSingular );
    #include "El/core/NestedGuardAndPayload.h"
}

#define DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<T,W,Z>& X, \
    bool checkIfSingular );

#define DIST_PROTO_INNER_REAL(T,U,V,W,Z) \
  DIST_PROTO_INNER_BASE(T,U,V,W,Z) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<T,U,V>& d, DistMatrix<Complex<T>,W,Z>& X, \
    bool checkIfSingular );

#define DIST_PROTO_BASE(T,U,V) \
  DIST_PROTO_INNER_BASE(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER_BASE(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER_BASE(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER_BASE(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER_BASE(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER_BASE(T,U,V,VR,  STAR);

#define DIST_PROTO_REAL(T,U,V) \
  DIST_PROTO_INNER_REAL(T,U,V,CIRC,CIRC); \
  DIST_PROTO_INNER_REAL(T,U,V,MC,  MR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,MC,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,MD,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,MR,  MC  ); \
  DIST_PROTO_INNER_REAL(T,U,V,MR,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,MD  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,MR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,VC  ); \
  DIST_PROTO_INNER_REAL(T,U,V,STAR,VR  ); \
  DIST_PROTO_INNER_REAL(T,U,V,VC,  STAR); \
  DIST_PROTO_INNER_REAL(T,U,V,VR,  STAR);

#define PROTO_BASE(T) \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<T>& d, Matrix<T>& X, bool checkIfSingular ); \
  template void DiagonalSolve \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<T>& d, AbstractDistMatrix<T>& X, \
    bool checkIfSingular ); \
  DIST_PROTO_BASE(T,CIRC,CIRC); \
  DIST_PROTO_BASE(T,MC,  MR  ); \
  DIST_PROTO_BASE(T,MC,  STAR); \
  DIST_PROTO_BASE(T,MD,  STAR); \
  DIST_PROTO_BASE(T,MR,  MC  ); \
  DIST_PROTO_BASE(T,MR,  STAR); \
  DIST_PROTO_BASE(T,STAR,MC  ); \
  DIST_PROTO_BASE(T,STAR,MD  ); \
  DIST_PROTO_BASE(T,STAR,MR  ); \
  DIST_PROTO_BASE(T,STAR,STAR); \
  DIST_PROTO_BASE(T,STAR,VC  ); \
  DIST_PROTO_BASE(T,STAR,VR  ); \
  DIST_PROTO_BASE(T,VC  ,STAR); \
  DIST_PROTO_BASE(T,VR  ,STAR);

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

PROTO_REAL(float);
PROTO_REAL(double);
PROTO_BASE(Complex<float>);
PROTO_BASE(Complex<double>);

} // namespace El
