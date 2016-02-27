/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace triang_eig {

/* Determine machine dependent parameters to control overflow
 *   Note: LAPACK uses more complicated parameters to handle 
 *   issues that can happen on Cray machines.
 */
template<typename Real>
inline pair<Real,Real>
OverflowParameters()
{
    const Real unfl = lapack::MachineSafeMin<Real>();
    const Real ovfl = lapack::MachineOverflowThreshold<Real>();
    const Real ulp  = lapack::MachinePrecision<Real>();
    Real smallNum = Max( unfl/ulp, Real(1)/(ovfl*ulp) );
    Real bigNum = Real(1)/smallNum;
    return pair<Real,Real>(smallNum,bigNum);
}

template<typename F>
void
DiagonalBlockSolve
(       Matrix<F>& U,
  const Matrix<F>& shifts,
        Matrix<F>& X,
        Matrix<F>& scales )
{
    typedef Base<F> Real;
  
    DEBUG_ONLY(
      CSE cse("triang_eig::DiagonalBlockSolve");
      if( shifts.Height() != X.Width() )
          LogicError("Incompatible number of shifts");
    )
    auto diag = GetDiagonal(U);
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();

    auto params = OverflowParameters<Real>();
    const Real smallNum = params.first;
    const Real bigNum = params.second;

    const Real oneHalf = Real(1)/Real(2);
    const Real oneQuarter = Real(1)/Real(4);
    
    // Default scale is 1
    Ones( scales, numShifts, 1 );

    // Compute infinity norms of columns of U (excluding diagonal)
    // TODO: scale cnorm if an entry is bigger than bigNum
    Matrix<Real> cnorm( n, 1 );
    cnorm.Set( 0, 0, Real(0) );
    for( Int j=1; j<n; ++j )
    {
        cnorm.Set( j, 0, MaxNorm( U(IR(0,j),IR(j)) ) );
    }

    // Iterate through RHS's
    for( Int j=1; j<numShifts; ++j )
    {
        // Initialize triangular system
        ShiftDiagonal( U, -shifts.Get(j,0) );
        auto xj = X( IR(0,Min(n,j)), IR(j) );

        // Determine largest entry of RHS
        Real xjMax = MaxNorm( xj );
        if( xjMax >= bigNum )
        {
            const Real s = oneHalf*bigNum/xjMax;
            xj *= s;
            xjMax *= s;
            scales.Set( j, 0, s*scales.Get(j,0) );
        }
        xjMax = Max( xjMax, 2*smallNum );

        // Estimate growth of entries in triangular solve
        //   Note: See "Robust Triangular Solves for Use in Condition
        //   Estimation" by Edward Anderson for explanation of bounds.
        Real invGi = 1/xjMax;
        Real invMi = invGi;
        for( Int i=Min(n,j)-1; i>=0; --i )
        {
            const Real absUii = SafeAbs( U.Get(i,i) );
            if( invGi<=smallNum || invMi<=smallNum || absUii<=smallNum )
            {
                invGi = 0;
                break;
            }
            invMi = Min( invMi, absUii*invGi );
            if( i > 0 )
            {
                invGi *= absUii/(absUii+cnorm.Get(i,0));
            }
        }
        invGi = Min( invGi, invMi );

        // Call TRSV if estimated growth is not too large
        if( invGi > smallNum )
        {
            blas::Trsv
            ( 'U', 'N', 'N', Min(n,j),
              U.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
        }
        // Perform backward substitution if estimated growth is large
        else
        {
            for( Int i=Min(n,j)-1; i>=0; --i )
            {
                // Perform division and check for overflow
                const Real absUii = SafeAbs( U.Get(i,i) );
                Real absXij = SafeAbs( xj.Get(i,0) );
                if( absUii > smallNum )
                {
                    if( absUii<=1 && absXij>=absUii*bigNum )
                    {
                        // Set overflowing entry to 0.5/U[i,i]
                        const Real s = oneHalf/absXij;
                        scales.Set( j, 0, s*scales.Get(j,0) );
                        xj *= s;
                        xjMax *= s;
                    }
                    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
                }
                else if( absUii > 0 )
                {
                    if( absXij >= absUii*bigNum )
                    {
                        // Set overflowing entry to bigNum/2
                        const Real s = oneHalf*absUii*bigNum/absXij;
                        scales.Set( j, 0, s*scales.Get(j,0) );
                        xj *= s;
                        xjMax *= s;
                    }
                    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
                }
                else
                {
                    // TODO: maybe this tolerance should be loosened to
                    //   | Xij | >= || A || * eps
                    if( absXij >= smallNum )
                    {
                        scales.Set( j, 0, F(0) );
                        Zero( xj );
                        xjMax = 0;
                    }
                    xj.Set( i, 0, F(1) );
                }

                if( i > 0 )
                {
                    // Check for possible overflows in AXPY
                    // Note: G(i+1) <= G(i) + | Xij | * cnorm(i)
                    absXij = SafeAbs( xj.Get(i,0) );
                    if( absXij >= 1 &&
                        cnorm.Get(i,0) >= (bigNum-xjMax)/absXij )
                    {
                        const Real s = oneQuarter/absXij;
                        scales.Set( j, 0, s*scales.Get(j,0) );
                        xj *= s;
                        xjMax *= s;
                        absXij *= s;
                    }
                    else if( absXij < 1 &&
                             absXij*cnorm.Get(i,0) >= bigNum-xjMax )
                    {
                        const Real s = oneQuarter;
                        scales.Set( j, 0, s*scales.Get(j,0) );
                        xj *= s;
                        xjMax *= s;
                        absXij *= s;
                    }
                    xjMax += absXij * cnorm.Get(i,0);

                    // AXPY
                    auto U01 = U( IR(0,i), IR(i) );
                    auto X1  = X( IR(0,i), IR(j) );
                    Axpy( -xj.Get(i,0), U01, X1 );
                }
            }
        }

        // Reset matrix diagonal
        SetDiagonal( U, diag );
    }
}
  
} // namespace triang_eig

template<typename F>
void TriangEig( Matrix<F>& U, Matrix<F>& X ) 
{
  
    DEBUG_ONLY(CSE cse("TriangEig"))
    const Int m = U.Height();

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    // Solve multi-shift triangular system
    Matrix<F> shifts, scales;
    GetDiagonal( U, shifts );
    SafeMultiShiftTrsm( LEFT, UPPER, NORMAL, F(1), U, shifts, X, scales );
    SetDiagonal( X, scales );

    // Normalize eigenvectors
    for( Int j=0; j<m; ++j )
    {
        auto xj = X( IR(0,j+1), IR(j) );
        Scale( 1/Nrm2(xj), xj );
    }
}
  
template<typename F>
void TriangEig
( const ElementalMatrix<F>& UPre, 
        ElementalMatrix<F>& XPre ) 
{
    DEBUG_ONLY(CSE cse("TriangEig"))
    const Int m = UPre.Height();
      
    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    // Solve multi-shift triangular system
    const Grid& g = U.Grid();
    DistMatrix<F,VR,STAR> shifts(g), scales(g);
    GetDiagonal( U, shifts );
    SafeMultiShiftTrsm( LEFT, UPPER, NORMAL, F(1), U, shifts, X, scales );
    SetDiagonal( X, scales );
    
    // Normalize eigenvectors
    for( Int j=1; j<m; ++j )
    {
        auto xj = X( IR(0,j+1), IR(j) );
        xj *= 1/Nrm2(xj);
    }
}

#define PROTO(F) \
  template void TriangEig \
  (       Matrix<F>& T, \
          Matrix<F>& X ); \
  template void TriangEig \
  ( const ElementalMatrix<F>& T, \
          ElementalMatrix<F>& X );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
