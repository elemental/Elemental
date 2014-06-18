/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega )
{
    DEBUG_ONLY(CallStackEntry cse("FoxLi"))
    typedef Complex<Real> C;
    const Real pi = 4*Atan( Real(1) );
    const C phi = Sqrt( C(0,omega/pi) ); 
    
    // Compute Gauss quadrature points and weights
    Matrix<Real> d, e; 
    Zeros( d, n, 1 );
    e.Resize( n-1, 1 );
    for( Int j=0; j<n-1; ++j )
    {
        const Real betaInv = 2*Sqrt(1-Pow(j+Real(1),-2)/4);
        e.Set( j, 0, 1/betaInv );
    }
    Matrix<Real> x, Z;
    HermitianTridiagEig( d, e, x, Z, UNSORTED );
    auto z = LockedView( Z, 0, 0, 1, n );
    Matrix<Real> sqrtWeights( z ), sqrtWeightsTrans;
    for( Int j=0; j<n; ++j )
        sqrtWeights.Set( 0, j, Sqrt(2)*Abs(sqrtWeights.Get(0,j)) );
    herm_eig::Sort( x, sqrtWeights, ASCENDING );
    Transpose( sqrtWeights, sqrtWeightsTrans );

    // Form the integral operator
    A.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            const Real theta = -omega*Pow(x.Get(i,0)-x.Get(j,0),2);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.Set( i, j, phi*C(realPart,imagPart) );
        }
    }

    // Apply the weighting
    DiagonalScale( LEFT, NORMAL, sqrtWeightsTrans, A );
    DiagonalScale( RIGHT, NORMAL, sqrtWeightsTrans, A );
}

template<typename Real,Dist U,Dist V>
void FoxLi( DistMatrix<Complex<Real>,U,V>& A, Int n, Real omega )
{
    DEBUG_ONLY(CallStackEntry cse("FoxLi"))
    typedef Complex<Real> C;
    const Real pi = 4*Atan( Real(1) );
    const C phi = Sqrt( C(0,omega/pi) ); 
    
    // Compute Gauss quadrature points and weights
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> d(g), e(g); 
    Zeros( d, n, 1 );
    e.Resize( n-1, 1 );
    for( Int iLoc=0; iLoc<e.LocalHeight(); ++iLoc )
    {
        const Int i = e.GlobalRow(iLoc);
        const Real betaInv = 2*Sqrt(1-Pow(i+Real(1),-2)/4);
        e.SetLocal( iLoc, 0, 1/betaInv );
    }
    DistMatrix<Real,VR,STAR> x(g);
    DistMatrix<Real,STAR,VR> Z(g);
    HermitianTridiagEig( d, e, x, Z, UNSORTED );
    auto z = LockedView( Z, 0, 0, 1, n );
    DistMatrix<Real,STAR,VR> sqrtWeights( z );
    for( Int jLoc=0; jLoc<sqrtWeights.LocalWidth(); ++jLoc )
        sqrtWeights.SetLocal
        ( 0, jLoc, Sqrt(2)*Abs(sqrtWeights.GetLocal(0,jLoc)) );
    herm_eig::Sort( x, sqrtWeights, ASCENDING );

    // Form the integral operator
    A.Resize( n, n );
    const Dist UGath = GatheredDist<U>();
    const Dist VGath = GatheredDist<V>();
    DistMatrix<Real,U,VGath> x_U_VGath( A.Grid() );
    DistMatrix<Real,V,UGath> x_V_UGath( A.Grid() );
    x_U_VGath.AlignWith( A ); 
    x_V_UGath.AlignWith( A );
    x_U_VGath = x;
    x_V_UGath = x;
    for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
    {
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
        {
            const Real diff = x_U_VGath.GetLocal(iLoc,0)-
                              x_V_UGath.GetLocal(jLoc,0);
            const Real theta = -omega*Pow(diff,2);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A.SetLocal( iLoc, jLoc, phi*C(realPart,imagPart) );
        }
    }

    // Apply the weighting
    DistMatrix<Real,VR,STAR> sqrtWeightsTrans(g);
    Transpose( sqrtWeights, sqrtWeightsTrans );
    DiagonalScale( LEFT, NORMAL, sqrtWeightsTrans, A );
    DiagonalScale( RIGHT, NORMAL, sqrtWeightsTrans, A );
}

#define PROTO_DIST(Real,U,V) \
  template void FoxLi( DistMatrix<Complex<Real>,U,V>& A, Int n, Real omega ); 

#define PROTO(Real) \
  template void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega ); \
  PROTO_DIST(Real,CIRC,CIRC) \
  PROTO_DIST(Real,MC,  MR  ) \
  PROTO_DIST(Real,MC,  STAR) \
  PROTO_DIST(Real,MD,  STAR) \
  PROTO_DIST(Real,MR,  MC  ) \
  PROTO_DIST(Real,MR,  STAR) \
  PROTO_DIST(Real,STAR,MC  ) \
  PROTO_DIST(Real,STAR,MD  ) \
  PROTO_DIST(Real,STAR,MR  ) \
  PROTO_DIST(Real,STAR,STAR) \
  PROTO_DIST(Real,STAR,VC  ) \
  PROTO_DIST(Real,STAR,VR  ) \
  PROTO_DIST(Real,VC,  STAR) \
  PROTO_DIST(Real,VR,  STAR)

//PROTO(float)
PROTO(double)

} // namespace El
