/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// This implementation is based upon the paper:
//
// Ivan Morel, Damien Stehle, and Gilles Villard,
// "H-LLL: Using Householder inside LLL",
// ISSAC '09, July 28--31, 2009, Seoul, Republic of Korea.

namespace El {

template<typename Real>
inline void Round( Matrix<Real>& B )
{
    DEBUG_ONLY(CSE cse("Round"))
    const Int m = B.Height();
    const Int n = B.Width();
    Real* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            BBuf[i+j*BLDim] = Round(BBuf[i+j*BLDim]);
}

namespace lll {

// NOTE: It is assumed that {delta,eta,theta} are valid
template<typename Real>
void Incomplete
( Int k,
  Matrix<Real>& B,
  Matrix<Real>& H,
  Matrix<Real>& t,
  Matrix<Real>& d,
  Real delta,
  Real eta,
  Real theta,
  Real loopTol,
  Real zeroTol )
{
    DEBUG_ONLY(CSE cse("lll::Incomplete"))
    const Int m = B.Height();

    vector<Real> xBuf(k);

    Real* BBuf = B.Buffer();
    const Int BLDim = B.LDim();

    Real* HBuf = H.Buffer();
    const Int HLDim = H.LDim();

    Real* dBuf = d.Buffer();
    Real* tBuf = t.Buffer();

    const bool progress = true;
    const bool print = false;

    while( true ) 
    {
        if( print )
        {
            Print( B, "B" );
            Print( H, "H" );
        }
        // Compute the (scaled) top of r_k from b_k using Householder
        // TODO: Maintain the reflectors in an accumulated form
        for( Int i=0; i<m; ++i )
            HBuf[i+k*HLDim] = BBuf[i+k*BLDim];
        for( Int i=0; i<k; ++i )
        {
            // Apply the i'th Householder reflector

            // Temporarily replace H(i,i) with 1
            const Real alpha = HBuf[i+i*HLDim]; 
            HBuf[i+i*HLDim] = 1;

            const Real innerProd =
              blas::Dot( m-i, &HBuf[i+i*HLDim], 1, &HBuf[i+k*HLDim], 1 );
            blas::Axpy
            ( m-i, -tBuf[i]*innerProd,
              &HBuf[i+i*HLDim], 1, &HBuf[i+k*HLDim], 1 );

            // Restore H(i,i)
            HBuf[i+i*HLDim] = alpha; 
        }
        // Fix the scaling of r_k
        for( Int i=0; i<k; ++i )
            HBuf[i+k*HLDim] *= dBuf[i];

        for( Int i=k-1; i>=0; --i )
        {
            const Real chi = Round(HBuf[i+k*HLDim]/HBuf[i+i*HLDim]);
            xBuf[i] = chi;
            blas::Axpy( i, -chi, &HBuf[i*HLDim], 1, &HBuf[k*HLDim], 1 );
        }

        const Real oldNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        blas::Gemv
        ( 'N', m, k,
          Real(-1), BBuf, BLDim, &xBuf[0], 1,
          Real(+1), &BBuf[k*BLDim], 1 );
        const Real newNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( newNorm*newNorm > loopTol*oldNorm*oldNorm )
        {
            break;
        }
        else if( progress )
            Output
            ("  Reorthogonalizing with k=",k,
             " since oldNorm=",oldNorm," and newNorm=",newNorm);
    }

    // Recompute {r_k,v_k,tau_k,delta_k} from b_k using Householder
    // TODO: Maintain the reflectors in an accumulated form
    for( Int i=0; i<m; ++i )
        HBuf[i+k*HLDim] = BBuf[i+k*BLDim];
    for( Int i=0; i<k; ++i )
    {
        // Apply the i'th Householder reflector

        // Temporarily replace H(i,i) with 1
        const Real alpha = HBuf[i+i*HLDim]; 
        HBuf[i+i*HLDim] = 1;

        const Real innerProd =
          blas::Dot( m-i, &HBuf[i+i*HLDim], 1, &HBuf[i+k*HLDim], 1 );
        blas::Axpy
        ( m-i, -tBuf[i]*innerProd, &HBuf[i+i*HLDim], 1, &HBuf[i+k*HLDim], 1 );

        // Restore H(i,i)
        HBuf[i+i*HLDim] = alpha; 
    }
    // Fix the scaling of r_k
    for( Int i=0; i<k; ++i )
        HBuf[i+k*HLDim] *= dBuf[i];
    // Perform the next step of Householder reduction
    {
        Real& rhokk = HBuf[k+k*HLDim]; 
        auto a21 = H( IR(k+1,END), IR(k) );
        tBuf[k] = LeftReflector( rhokk, a21 );
        if( rhokk < Real(0) )
        {
            dBuf[k] = -1;
            rhokk *= -1;
        }
        else
            dBuf[k] = 1;
        if( rhokk < zeroTol )
            throw SingularMatrixException();
    }
}

} // namespace lll

template<typename Real>
void LLL( Matrix<Real>& B, Real delta, Real eta, Real theta, Real loopTol )
{
    DEBUG_ONLY(CSE cse("LLL"))
    if( delta > Real(1) )
        LogicError("delta is assumed to be at most 1");
    if( eta < Real(1)/Real(2) )
        LogicError("eta is assumed to be at least 1/2");
    if( theta >= eta - Real(1)/Real(2) )
        LogicError("We assume that theta < eta - 1/2");

    const bool progress = true;
    const bool print = false;
    Real zeroTol = Pow(Epsilon<Real>(),Real(0.5));

    // Force the input to be integer-valued; it would be okay to assume this
    Round( B );
    if( print )
        Print( B, "Round(B)" );

    const Int m = B.Height();
    const Int n = B.Width();
    const Int minDim = Min(m,n);
    Matrix<Real> H, t, d;
    Zeros( H, m, n );
    Zeros( d, minDim, 1 );
    Zeros( t, minDim, 1 );

    // Perform the first step of Householder reduction
    {
        Real* HBuf = H.Buffer();
        Real* BBuf = B.Buffer();
        for( Int i=0; i<m; ++i )
            HBuf[i] = BBuf[i];
        auto alpha11 = H( IR(0), IR(0) );
        auto a21 = H( IR(1,END), IR(0) );
        const Real tau = LeftReflector( alpha11, a21 );
        t.Set( 0, 0, tau );
        if( HBuf[0] < Real(0) )
        {
            d.Set(0,0,Real(-1));
            HBuf[0] *= Real(-1);
        }
        else
            d.Set(0,0,Real(1));
        if( HBuf[0] < zeroTol )
            throw SingularMatrixException();
    }

    Int k=1;
    while( k < n )
    {
        lll::Incomplete( k, B, H, t, d, delta, eta, theta, loopTol, zeroTol );
        const Real bNorm = FrobeniusNorm( B(ALL,IR(k)) );
        const Real rTNorm = FrobeniusNorm( H(IR(0,k-1),IR(k)) );
        const Real s = bNorm*bNorm - rTNorm*rTNorm;
        const Real rho_km1_km1 = H.Get(k-1,k-1);
        if( delta*rho_km1_km1*rho_km1_km1 <= s )
        {
            ++k;
        }
        else
        {
            if( progress )
                Output("Dropping from k=",k," to ",Max(k-1,1));
            ColSwap( B, k-1, k );
            k = Max(k-1,1);
        }
    }
}

#define PROTO(Real) \
  template void LLL \
  ( Matrix<Real>& B, Real delta, Real eta, Real theta, Real loopTol );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
