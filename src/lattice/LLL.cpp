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
//
// Note that there is a subtle bug in Algorithm 4 of the paper, as step 7
// involves setting k := max(k-1,2), which, in the case of zero indexing,
// becomes k := max(k-1,1). The second branch of the 'max' is only activated
// in the case where k=1 and is kept at this value despite having just 
// swapped b_0 and b_1. But this would imply that the 0'th column of the 
// Householder QR factorization is now incorrect, and so it must be refreshed.
// This implementation fixes said issue via a call to lll::HouseholderStep
// with k=0.
//
// Future work will involve investigating blocked algorithms and/or 
// distributed-memory and/or GPU implementations.
//
// Also, the following paper was consulted for more general factorization
// viewpoints:
//
// Dirk Wubben, Ronald Bohnke, Volker Kuhn, and Karl-Dirk Kammeyer,
// "MMSE-Based Lattice-Reduction for Near-ML Detection of MIMO Systems",
// ITG Workshop on Smart Antennas, pp. 106--113, 2004

namespace El {

namespace lll {

// Put the k'th column of B into the k'th column of QR and then rotate
// said column with the first k-1 (scaled) Householder reflectors.
//
// TODO: Maintain the reflectors in an accumulated form
template<typename F>
void ExpandQR
( Int k,
  const Matrix<F>& B,
        Matrix<F>& QR,
  const Matrix<F>& t,
  const Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CSE cse("lll::ExpandQR"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const F* BBuf = B.LockedBuffer();
          F* QRBuf = QR.Buffer();  
    const Int BLDim = B.LDim();
    const Int QRLDim = QR.LDim();

    // Copy in the k'th column of B
    for( Int i=0; i<m; ++i )
        QRBuf[i+k*QRLDim] = BBuf[i+k*BLDim];

    for( Int i=0; i<k; ++i )
    {
        // Apply the i'th Householder reflector

        // Temporarily replace QR(i,i) with 1
        const Real alpha = RealPart(QRBuf[i+i*QRLDim]);
        QRBuf[i+i*QRLDim] = 1;

        const F innerProd =
          blas::Dot
          ( m-i,
            &QRBuf[i+i*QRLDim], 1,
            &QRBuf[i+k*QRLDim], 1 );
        blas::Axpy
        ( m-i, -t.Get(i,0)*innerProd,
          &QRBuf[i+i*QRLDim], 1,
          &QRBuf[i+k*QRLDim], 1 );

        // Fix the scaling
        QRBuf[i+k*QRLDim] *= d.Get(i,0);

        // Restore H(i,i)
        QRBuf[i+i*QRLDim] = alpha; 
    }
}

template<typename F>
void HouseholderStep
( Int k,
  const Matrix<F>& B,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
        Base<F> zeroTol )
{
    DEBUG_ONLY(CSE cse("lll::HouseholderStep"))
    typedef Base<F> Real;

    lll::ExpandQR( k, B, QR, t, d );

    F* QRBuf = QR.Buffer();
    const Int QRLDim = QR.LDim();

    // Perform the next step of Householder reduction
    F& rhokk = QRBuf[k+k*QRLDim]; 
    auto qr21 = QR( IR(k+1,END), IR(k) );
    F tau = LeftReflector( rhokk, qr21 );
    t.Set( k, 0, tau );
    if( RealPart(rhokk) < Real(0) )
    {
        d.Set( k, 0, -1 );
        rhokk *= -1;
    }
    else
        d.Set( k, 0, +1 );
    if( RealPart(rhokk) < zeroTol )
        throw SingularMatrixException();
}

// NOTE: It is assumed that delta is valid
template<typename F>
void Step
( Int k,
  Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  Base<F> delta,
  Base<F> loopTol,
  Base<F> zeroTol,
  bool progress )
{
    DEBUG_ONLY(CSE cse("lll::Step"))
    typedef Base<F> Real;
    const Int m = B.Height();

    vector<F> xBuf(k);

    F* BBuf = B.Buffer();
    F* QRBuf = QR.Buffer();
    const Int BLDim = B.LDim();
    const Int QRLDim = QR.LDim();

    while( true ) 
    {
        lll::ExpandQR( k, B, QR, t, d );

        for( Int i=k-1; i>=0; --i )
        {
            const F chi = Round(QRBuf[i+k*QRLDim]/QRBuf[i+i*QRLDim]);
            xBuf[i] = chi;
            blas::Axpy( i, -chi, &QRBuf[i*QRLDim], 1, &QRBuf[k*QRLDim], 1 );
        }

        const Real oldNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        blas::Gemv
        ( 'N', m, k,
          F(-1), BBuf, BLDim, &xBuf[0], 1,
          F(+1), &BBuf[k*BLDim], 1 );
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

    lll::HouseholderStep( k, B, QR, t, d, zeroTol );
}

} // namespace lll

template<typename F>
Int LLL
( Matrix<F>& B,
  Matrix<F>& QR,
  Base<F> delta,
  Base<F> loopTol,
  bool presort,
  bool smallestFirst,
  bool progress )
{
    DEBUG_ONLY(CSE cse("LLL"))
    typedef Base<F> Real;
    if( delta > Real(1) )
        LogicError("delta is assumed to be at most 1");

    Real zeroTol = Pow(Epsilon<Real>(),Real(0.5));

    // Force the input to be integer-valued; it would be okay to assume this
    Round( B );

    if( presort )
    {
        QRCtrl<Real> ctrl;
        ctrl.smallestFirst = smallestFirst;

        auto BCopy = B;
        Matrix<F> t;
        Matrix<Real> d;
        Matrix<Int> colPerm;
        El::QR( BCopy, t, d, colPerm, ctrl );

        InversePermuteCols( B, colPerm );
    }

    const Int m = B.Height();
    const Int n = B.Width();
    const Int minDim = Min(m,n);
    Matrix<F> t;
    Matrix<Real> d;
    Zeros( QR, m, n );
    Zeros( d, minDim, 1 );
    Zeros( t, minDim, 1 );

    // Perform the first step of Householder reduction
    lll::HouseholderStep( 0, B, QR, t, d, zeroTol );

    Int k=1, numBacktrack=0;
    while( k < n )
    {
        lll::Step( k, B, QR, t, d, delta, loopTol, zeroTol, progress );

        const Real bNorm = FrobeniusNorm( B(ALL,IR(k)) );
        const Real rTNorm = FrobeniusNorm( QR(IR(0,k-1),IR(k)) );
        const Real s = bNorm*bNorm - rTNorm*rTNorm;
        const Real rho_km1_km1 = QR.GetRealPart(k-1,k-1);
        if( delta*rho_km1_km1*rho_km1_km1 <= s )
        {
            ++k;
        }
        else
        {
            ++numBacktrack;
            if( progress )
                Output("Dropping from k=",k," to ",Max(k-1,1));
            ColSwap( B, k-1, k );
            if( k == 1 )
            {
                // We must reinitialize since we keep k=1
                lll::HouseholderStep( 0, B, QR, t, d, zeroTol );
            }
            else
            {
                k = k-1; 
            }
        }
    }
    return numBacktrack;
}

template<typename F>
Base<F> LLLDelta( const Matrix<F>& QR )
{
    DEBUG_ONLY(CSE cse("LLLDelta"))
    typedef Base<F> Real;
    const Int m = QR.Height();
    const Int n = QR.Width();
    const Int minDim = Min(m,n);

    auto QRTop = QR( IR(0,minDim), ALL );
    auto R = QRTop;
    MakeTrapezoidal( UPPER, R );
    
    // Find the maximum delta such that
    //
    //   delta R(k,k)^2 <= R(k+1,k+1)^2 + |R(k,k+1)|^2
    //
    // for 0 <= k < n-1.
    //
    if( n <= 1 )
        return 1; // the best-possible delta
    Matrix<F> z;
    Real delta = std::numeric_limits<Real>::max();
    for( Int i=0; i<n-1; ++i )
    {
        const Real rho_i_i = R.GetRealPart(i,i);
        const Real rho_i_ip1 = Abs(R.Get(i,i+1));
        const Real rho_ip1_ip1 = R.GetRealPart(i+1,i+1);

        const Real deltaBound =
          (rho_ip1_ip1*rho_ip1_ip1+rho_i_ip1*rho_i_ip1)/(rho_i_i*rho_i_i);
 
        delta = Min(delta,deltaBound);
    }

    // Ensure that
    //
    //    | R(l,k) | <= 0.5 | R(l,l) | for all 0 <= l < k < n
    //
    // NOTE: This does not seem to hold for complex LLL reductions,
    //       so sqrt(2)/2 is used instead.
    auto diagR = GetDiagonal(R);
    DiagonalSolve( LEFT, NORMAL, diagR, R );
    ShiftDiagonal( R, F(-1) );
    const Real maxRatio = MaxNorm( R );
    const Real bound = 
      ( IsComplex<F>::value ? Real(1)/Sqrt(Real(2)) : Real(1)/Real(2) );
    if( maxRatio > bound+Pow(Epsilon<Real>(),Real(3)/Real(4)) )
        return 0; // the worst-possible delta

    return delta;
}

#define PROTO(F) \
  template Int LLL \
  ( Matrix<F>& B, \
    Matrix<F>& QR, \
    Base<F> delta, \
    Base<F> loopTol, \
    bool presort, \
    bool smallestFirst, \
    bool progress ); \
  template Base<F> LLLDelta( const Matrix<F>& QR );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
