/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("FoxLi"))
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
    auto z = Z( IR(0), ALL );
    Matrix<Real> sqrtWeights( z ), sqrtWeightsTrans;
    for( Int j=0; j<n; ++j )
        sqrtWeights.Set( 0, j, Sqrt(Real(2))*Abs(sqrtWeights.Get(0,j)) );
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

template<typename Real>
void FoxLi( ElementalMatrix<Complex<Real>>& APre, Int n, Real omega )
{
    DEBUG_ONLY(CSE cse("FoxLi"))
    typedef Complex<Real> C;
    const Real pi = 4*Atan( Real(1) );
    const C phi = Sqrt( C(0,omega/pi) ); 

    DistMatrixWriteProxy<C,C,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    
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
    auto z = Z( IR(0), ALL );
    DistMatrix<Real,STAR,VR> sqrtWeights( z );
    for( Int jLoc=0; jLoc<sqrtWeights.LocalWidth(); ++jLoc )
        sqrtWeights.SetLocal
            ( 0, jLoc, Sqrt(Real(2))*Abs(sqrtWeights.GetLocal(0,jLoc)) );
    herm_eig::Sort( x, sqrtWeights, ASCENDING );

    // Form the integral operator
    A.Resize( n, n );
    DistMatrix<Real,MC,STAR> x_MC_STAR( A.Grid() );
    DistMatrix<Real,MR,STAR> x_MR_STAR( A.Grid() );
    x_MC_STAR.AlignWith( A ); 
    x_MR_STAR.AlignWith( A );
    x_MC_STAR = x;
    x_MR_STAR = x;
    for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
    {
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
        {
            const Real diff = x_MC_STAR.GetLocal(iLoc,0)-
                              x_MR_STAR.GetLocal(jLoc,0);
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

#define PROTO(Real) \
  template void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega ); \
  template void FoxLi \
  ( ElementalMatrix<Complex<Real>>& A, Int n, Real omega );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
