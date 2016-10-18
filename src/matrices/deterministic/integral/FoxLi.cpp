/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Real>
void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega )
{
    DEBUG_CSE
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
        e(j) = 1/betaInv;
    }
    Matrix<Real> x, Z;
    HermitianTridiagEigCtrl<Real> ctrl;
    ctrl.sort = UNSORTED;
    HermitianTridiagEig( d, e, x, Z, ctrl );
    auto z = Z( IR(0), ALL );
    Matrix<Real> sqrtWeights( z ), sqrtWeightsTrans;
    for( Int j=0; j<n; ++j )
        sqrtWeights(0,j) = Sqrt(Real(2))*Abs(sqrtWeights(0,j));
    auto sortPairs = TaggedSort( x, ASCENDING );
    for( Int j=0; j<n; ++j )
        x(j) = sortPairs[j].value;
    ApplyTaggedSortToEachRow( sortPairs, sqrtWeights );
    Transpose( sqrtWeights, sqrtWeightsTrans );

    // Form the integral operator
    A.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            const Real theta = -omega*Pow(x(i)-x(j),2);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            A(i,j) = phi*C(realPart,imagPart);
        }
    }

    // Apply the weighting
    DiagonalScale( LEFT, NORMAL, sqrtWeightsTrans, A );
    DiagonalScale( RIGHT, NORMAL, sqrtWeightsTrans, A );
}

template<typename Real>
void FoxLi( ElementalMatrix<Complex<Real>>& APre, Int n, Real omega )
{
    DEBUG_CSE
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
    auto& eLoc = e.Matrix();
    for( Int iLoc=0; iLoc<e.LocalHeight(); ++iLoc )
    {
        const Int i = e.GlobalRow(iLoc);
        const Real betaInv = 2*Sqrt(1-Pow(i+Real(1),-2)/4);
        eLoc(iLoc) = 1/betaInv;
    }
    DistMatrix<Real,VR,STAR> x(g);
    DistMatrix<Real,STAR,VR> Z(g);
    HermitianTridiagEigCtrl<Real> ctrl;
    ctrl.sort = UNSORTED;
    HermitianTridiagEig( d, e, x, Z, ctrl );
    auto z = Z( IR(0), ALL );
    DistMatrix<Real,STAR,VR> sqrtWeights( z );
    auto& sqrtWeightsLoc = sqrtWeights.Matrix();
    for( Int jLoc=0; jLoc<sqrtWeights.LocalWidth(); ++jLoc )
        sqrtWeightsLoc(0,jLoc) = Sqrt(Real(2))*Abs(sqrtWeightsLoc(0,jLoc));
    auto sortPairs = TaggedSort( x, ASCENDING );
    for( Int j=0; j<n; ++j )
        x.Set( j, 0, sortPairs[j].value );
    ApplyTaggedSortToEachRow( sortPairs, sqrtWeights );

    // Form the integral operator
    A.Resize( n, n );
    DistMatrix<Real,MC,STAR> x_MC( A.Grid() );
    DistMatrix<Real,MR,STAR> x_MR( A.Grid() );
    x_MC.AlignWith( A ); 
    x_MR.AlignWith( A );
    x_MC = x;
    x_MR = x;
    auto& ALoc = A.Matrix();
    auto& x_MCLoc = x_MC.Matrix();
    auto& x_MRLoc = x_MR.Matrix();
    for( Int jLoc=0; jLoc<A.LocalWidth(); ++jLoc )
    {
        for( Int iLoc=0; iLoc<A.LocalHeight(); ++iLoc )
        {
            const Real diff = x_MCLoc(iLoc)-x_MRLoc(jLoc);
            const Real theta = -omega*Pow(diff,2);
            const Real realPart = Cos(theta);
            const Real imagPart = Sin(theta);
            ALoc(iLoc,jLoc) = phi*C(realPart,imagPart);
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
#define EL_ENABLE_QUAD
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
