/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// The implementations of Householder-based LLL in this file are extensions of
// the algorithm discussed in:
//
//   Ivan Morel, Damien Stehle, and Gilles Villard,
//   "H-LLL: Using Householder inside LLL",
//   ISSAC '09, July 28--31, 2009, Seoul, Republic of Korea.
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
// Furthermore, an analogue of the above algorithm which maintains an 
// accumulated set of Householder transformations, using a so-called 
// UT Transform,
//
//   Thierry Joffrain, Tze Meng Low, Enrique S. Quintana-Orti, and 
//   Robert van de Geijn,
//   "Accumulating Householder Transforms, Revisited",
//   ACM Transactions on Mathematical Software, Vol. 32, No. 2, pp. 169--179,
//   2006.
//
// However, it can be seen that the k'th iteration of the accumulated algorithm
// requires
//
//     6 k ( m - k/2 )
//
// operations, whereas the standard algorithm only requires
//
//     4 k ( m - k/2 )
// 
// operations. For this reason, there is not a clear benefit to using the
// accumulated algorithm for LLL in a sequential setting since the transition
// from level 1 to level 2 BLAS does not warrant the 50% increase in work.
// However, the situation might be significantly different in a distributed
// setting, where each inner product would involve a non-trivial amount of
// latency. Also, future work on 'greedy' variants of LLL may further swing
// the tide towards accumulated transforms.
//
// Also, the following paper was consulted for more general factorization
// viewpoints:
//
// Dirk Wubben, Ronald Bohnke, Volker Kuhn, and Karl-Dirk Kammeyer,
// "MMSE-Based Lattice-Reduction for Near-ML Detection of MIMO Systems",
// ITG Workshop on Smart Antennas, pp. 106--113, 2004
//
// Future work will involve investigating blocked algorithms and/or 
// distributed-memory and/or GPU implementations.
//
// The seminal work on distributed-memory implementations of LLL is
//
//   G. Villard, "Parallel lattice basis reduction",
//   Papers from the International Symposium on Symbolic and algebraic
//   computation, 1992, pp. 269--277.
//
// The key idea is that, because only |R(i,i+1)/R(i,i)| <= 1/2 is required
// in order to achieve the optimality bound of the first column of the 
// reduced basis (noticed by Lenstra et al. in "Factoring polynomials with
// rational coefficients"), one can transition to an "all swaps" algorithm
// which alternates between swapping all possible (b_i,b_{i+1}) columns with
// odd and even coefficients. The result is then only weakly LLL reduced,
// but the parallel algorithm is greatly simplified. Furthermore, the QR
// factorization can easily be patched up with an independent Givens rotation
// for each swap.
//
// Henri Cohen's "A course in computational algebraic number theory" was also
// heavily consulted for extending LLL to support linearly-dependent columns
// (and its extension to computing integer kernels).

namespace El {

static Timer applyHouseTimer, roundTimer, formSInvTimer;

namespace lll {

// Put the k'th column of B into the k'th column of QR and then rotate
// said column with the first k-1 (scaled) Householder reflectors.

template<typename F>
void ExpandQR
( Int k,
  const Matrix<F>& B,
        Matrix<F>& QR,
  const Matrix<F>& t,
  const Matrix<Base<F>>& d,
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::ExpandQR"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const F* BBuf = B.LockedBuffer();
          F* QRBuf = QR.Buffer();  
    const F* tBuf = t.LockedBuffer();
    const Base<F>* dBuf = d.LockedBuffer();
    const Int BLDim = B.LDim();
    const Int QRLDim = QR.LDim();

    // Copy in the k'th column of B
    for( Int i=0; i<m; ++i )
        QRBuf[i+k*QRLDim] = BBuf[i+k*BLDim];

    if( time )
        applyHouseTimer.Start();
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
        ( m-i, -tBuf[i]*innerProd,
          &QRBuf[i+i*QRLDim], 1,
          &QRBuf[i+k*QRLDim], 1 );

        // Fix the scaling
        QRBuf[i+k*QRLDim] *= dBuf[i];

        // Restore H(i,i)
        QRBuf[i+i*QRLDim] = alpha; 
    }
    if( time )
        applyHouseTimer.Stop();
}

template<typename F>
void HouseholderStep
( Int k,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::HouseholderStep"))
    typedef Base<F> Real;

    // Perform the next step of Householder reduction
    F* QRBuf = QR.Buffer();
    const Int QRLDim = QR.LDim();
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
}

// Return true if the new column is a zero vector
template<typename F>
bool Step
( Int k,
  Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool formU,
  bool formUInv,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::Step"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();

    vector<F> xBuf(k);

    F* BBuf = B.Buffer();
    F* UBuf = U.Buffer();
    F* UInvBuf = UInv.Buffer();
    F* QRBuf = QR.Buffer();
    const Int BLDim = B.LDim();
    const Int ULDim = U.LDim();
    const Int UInvLDim = UInv.LDim();
    const Int QRLDim = QR.LDim();

    while( true ) 
    {
        lll::ExpandQR( k, B, QR, t, d, ctrl.time );

        const Real oldNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( oldNorm <= ctrl.zeroTol )
        {
            for( Int i=0; i<m; ++i )
                BBuf[i+k*BLDim] = 0;
            for( Int i=0; i<m; ++i )
                QRBuf[i+k*QRLDim] = 0;
            t.Set( k, 0, Real(1)/Real(2) );
            d.Set( k, 0, Real(1) );
            return true;
        }

        if( ctrl.time )
            roundTimer.Start();
        if( ctrl.weak )
        {
            if( Abs(QRBuf[(k-1)+(k-1)*QRLDim]) > ctrl.zeroTol )
            {
                const F chi =
                  Round(QRBuf[(k-1)+k*QRLDim]/QRBuf[(k-1)+(k-1)*QRLDim]);
                xBuf[k-1] = chi;
                if( Abs(chi) > Real(1)/Real(2) )
                {
                    blas::Axpy
                    ( k, -chi,
                      &QRBuf[(k-1)*QRLDim], 1,
                      &QRBuf[k*QRLDim], 1 );
                    blas::Axpy
                    ( m, -chi,
                      &BBuf[(k-1)*BLDim], 1,
                      &BBuf[k*BLDim], 1 );
                    if( formU )
                        blas::Axpy
                        ( n, -chi,
                          &UBuf[(k-1)*ULDim], 1,
                          &UBuf[k*ULDim], 1 );
                    if( formUInv )
                        blas::Axpy
                        ( n, chi,
                          &UInvBuf[k], UInvLDim,
                          &UInvBuf[k-1], UInvLDim );
                }
            }
        }
        else
        {
            for( Int i=k-1; i>=0; --i )
            {
                if( Abs(QRBuf[i+i*QRLDim]) <= ctrl.zeroTol )
                {
                    xBuf[i] = 0;
                    continue;
                }
                const F chi = Round(QRBuf[i+k*QRLDim]/QRBuf[i+i*QRLDim]);
                xBuf[i] = chi;
                if( Abs(chi) > Real(1)/Real(2) )
                    blas::Axpy
                    ( i+1, -chi,
                      &QRBuf[i*QRLDim], 1,
                      &QRBuf[k*QRLDim], 1 );
            }
            blas::Gemv
            ( 'N', m, k,
              F(-1), BBuf, BLDim, &xBuf[0], 1,
              F(+1), &BBuf[k*BLDim], 1 );
            if( formU )
                blas::Gemv
                ( 'N', n, k,
                  F(-1), UBuf, ULDim, &xBuf[0], 1,
                  F(+1), &UBuf[k*ULDim], 1 );
            if( formUInv )
                blas::Geru
                ( k, n,
                  F(1), &xBuf[0],    1,
                        &UInvBuf[k], UInvLDim,
                        &UInvBuf[0], UInvLDim );
        }
        const Real newNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( ctrl.time )
            roundTimer.Stop();

        if( newNorm*newNorm > ctrl.reorthogTol*oldNorm*oldNorm )
        {
            break;
        }
        else if( ctrl.progress )
            Output
            ("  Reorthogonalizing with k=",k,
             " since oldNorm=",oldNorm," and newNorm=",newNorm);
    }

    lll::HouseholderStep( k, QR, t, d, ctrl.time );
    return false;
}

// Assume that V is m x n and SInv is n x n with, the first k columns of V
// and the first k rows and columns of SInv, up to date. 
template<typename F>
void ExpandBlockQR
( Int k,
  const Matrix<F>& B,
        Matrix<F>& QR,
        Matrix<F>& V,
        Matrix<F>& SInv,
  const Matrix<Base<F>>& d,
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::ExpandBlockQR"))
    const Int m = B.Height();
    const F* BBuf = B.LockedBuffer();
          F* QRBuf = QR.Buffer();
    const Base<F>* dBuf = d.LockedBuffer();
    const Int BLDim = B.LDim();
    const Int QRLDim = QR.LDim();

    // Copy in the k'th column of B
    for( Int i=0; i<m; ++i )
        QRBuf[i+k*QRLDim] = BBuf[i+k*BLDim];

    // Apply the first k Householder reflectors
    Matrix<F> z;

    if( time )
        applyHouseTimer.Start();
    // Exploit zeros in upper triangle of V
    /*
    Zeros( z, k, 1 );
    blas::Gemv
    ( 'C', m, k,
      F(1), V.Buffer(), V.LDim(),
            B.LockedBuffer(0,k), 1,
      F(1), z.Buffer(), 1 );
    */
    z.Resize( k, 1 );
    F* zBuf = z.Buffer();
    for( Int i=0; i<k; ++i )
        zBuf[i] = BBuf[i+k*BLDim];
    blas::Trmv( 'L', 'C', 'N', k, V.Buffer(), V.LDim(), z.Buffer(), 1 );
    blas::Gemv
    ( 'C', m-k, k,
      F(1), V.Buffer(k,0), V.LDim(),
            B.LockedBuffer(k,k), 1,
      F(1), z.Buffer(), 1 );

    blas::Trsv
    ( 'L', 'N', 'N', k, SInv.LockedBuffer(), SInv.LDim(), z.Buffer(), 1 );

    // Exploit zeros in upper triangle of V
    /*
    blas::Gemv
    ( 'N', m, k,
      F(-1), V.Buffer(), V.LDim(),
             z.LockedBuffer(), 1,
      F(1), QR.Buffer(0,k), 1 );
    */
    blas::Gemv
    ( 'N', m-k, k,
      F(-1), V.Buffer(k,0), V.LDim(),
             z.LockedBuffer(), 1,
      F(1), &QRBuf[k+k*QRLDim], 1 );
    blas::Trmv( 'L', 'N', 'N', k, V.Buffer(), V.LDim(), zBuf, 1 );
    for( Int i=0; i<k; ++i )
        QRBuf[i+k*QRLDim] -= zBuf[i];

    if( time )
        applyHouseTimer.Stop();

    // Fix the scaling
    for( Int i=0; i<k; ++i )
        QRBuf[i+k*QRLDim] *= dBuf[i];
}

// Thus, only V(:,k) and SInv(0:k,k) needs to be computed in order to 
// apply the block Householder transform I - V inv(S) V^H
template<typename F>
void BlockHouseholderStep
( Int k,
  Matrix<F>& QR,
  Matrix<F>& V,
  Matrix<F>& SInv,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::BlockHouseholderStep"))
    typedef Base<F> Real;
    const Int m = QR.Height();

    F* QRBuf = QR.Buffer();
    F* VBuf = V.Buffer();
    F* SInvBuf = SInv.Buffer();
    const Int QRLDim = QR.LDim();
    const Int VLDim = V.LDim();
    const Int SInvLDim = SInv.LDim();

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

    // Form the k'th column of V 
    for( Int i=0; i<k; ++i )
        VBuf[i+k*VLDim] = 0;
    VBuf[k+k*VLDim] = 1;
    for( Int i=k+1; i<m; ++i )
        VBuf[i+k*VLDim] = QRBuf[i+k*QRLDim];

    // Form the k'th row of SInv
    if( time )
        formSInvTimer.Start();
    /*
    blas::Gemv
    ( 'C', m, k,
      F(1), VBuf, VLDim, &VBuf[k*VLDim], 1,
      F(0), &SInvBuf[k], SInvLDim );
    */
    blas::Gemv
    ( 'C', m-k, k,
      F(1), &VBuf[k], VLDim, &VBuf[k+k*VLDim], 1,
      F(0), &SInvBuf[k], SInvLDim );
    SInvBuf[k+k*SInvLDim] = F(1)/t.Get(k,0);
    if( time )
        formSInvTimer.Stop();
}

// Return true if the new vector is a zero vector
template<typename F>
bool BlockStep
( Int k,
  Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  Matrix<F>& V,
  Matrix<F>& SInv,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool formU,
  bool formUInv,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::BlockStep"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();

    vector<F> xBuf(k);

    F* BBuf = B.Buffer();
    F* UBuf = U.Buffer();
    F* UInvBuf = UInv.Buffer();
    F* QRBuf = QR.Buffer();
    const Int BLDim = B.LDim();
    const Int ULDim = U.LDim();
    const Int UInvLDim = UInv.LDim();
    const Int QRLDim = QR.LDim();

    while( true ) 
    {
        lll::ExpandBlockQR( k, B, QR, V, SInv, d, ctrl.time );

        const Real oldNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( oldNorm <= ctrl.zeroTol )
        {
            for( Int i=0; i<m; ++i )
                BBuf[i+k*BLDim] = 0;
            for( Int i=0; i<m; ++i )
                QRBuf[i+k*QRLDim] = 0;
            t.Set( k, 0, Real(1)/Real(2) );
            d.Set( k, 0, Real(1) );
            return true;
        }

        if( ctrl.time )
            roundTimer.Start();
        if( ctrl.weak )
        {
            if( Abs(QRBuf[(k-1)+(k-1)*QRLDim]) > ctrl.zeroTol )
            {
                const F chi =
                  Round(QRBuf[(k-1)+k*QRLDim]/QRBuf[(k-1)+(k-1)*QRLDim]);
                xBuf[k-1] = chi;
                if( Abs(chi) > Real(1)/Real(2) )
                {
                    blas::Axpy
                    ( k, -chi,
                      &QRBuf[(k-1)*QRLDim], 1,
                      &QRBuf[k*QRLDim], 1 );
                    blas::Axpy
                    ( m, -chi,
                      &BBuf[(k-1)*BLDim], 1,
                      &BBuf[k*BLDim], 1 );
                    if( formU )
                        blas::Axpy
                        ( n, -chi,
                          &UBuf[(k-1)*ULDim], 1,
                          &UBuf[k*ULDim], 1 );
                    if( formUInv )
                        blas::Axpy
                        ( n, chi,
                          &UInvBuf[k], UInvLDim,
                          &UInvBuf[k-1], UInvLDim );
                }
            }
        }
        else
        {
            for( Int i=k-1; i>=0; --i )
            {
                if( Abs(QRBuf[i+i*QRLDim]) <= ctrl.zeroTol )
                {
                    xBuf[i] = 0;
                    continue;
                }
                const F chi = Round(QRBuf[i+k*QRLDim]/QRBuf[i+i*QRLDim]);
                xBuf[i] = chi;
                if( Abs(chi) > Real(1)/Real(2) )
                    blas::Axpy
                    ( i+1, -chi,
                      &QRBuf[i*QRLDim], 1,
                      &QRBuf[k*QRLDim], 1 );
            }
            blas::Gemv
            ( 'N', m, k,
              F(-1), BBuf, BLDim, &xBuf[0], 1,
              F(+1), &BBuf[k*BLDim], 1 );
            if( formU )
                blas::Gemv
                ( 'N', n, k,
                  F(-1), UBuf, ULDim, &xBuf[0], 1,
                  F(+1), &UBuf[k*ULDim], 1 );
            if( formUInv )
                blas::Geru
                ( k, n,
                  F(1), &xBuf[0],    1,
                        &UInvBuf[k], UInvLDim,
                        &UInvBuf[0], UInvLDim );
        }
        const Real newNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( ctrl.time )
            roundTimer.Stop();

        if( newNorm*newNorm > ctrl.reorthogTol*oldNorm*oldNorm )
        {
            break;
        }
        else if( ctrl.progress )
            Output
            ("  Reorthogonalizing with k=",k,
             " since oldNorm=",oldNorm," and newNorm=",newNorm);
    }

    lll::BlockHouseholderStep( k, QR, V, SInv, t, d, ctrl.time );
    return false;
}

template<typename F>
LLLInfo UnblockedAlg
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  bool formU,
  bool formUInv,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::UnblockedAlg"))
    typedef Base<F> Real;
    if( ctrl.time )
    {
        applyHouseTimer.Reset();
        roundTimer.Reset();
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
    lll::ExpandQR( 0, B, QR, t, d, ctrl.time );
    lll::HouseholderStep( 0, QR, t, d, ctrl.time );
    Int nullity = 0;
    {
        auto b0 = B(ALL,IR(0));
        if( FrobeniusNorm(b0) <= ctrl.zeroTol )
        {
            auto QR0 = QR(ALL,IR(0));
            Zero( b0 );
            Zero( QR0 );
            nullity = 1;
        }
    }

    Int k=1, numSwaps=0;
    while( k < n )
    {
        bool zeroVector =
          lll::Step( k, B, U, UInv, QR, t, d, formU, formUInv, ctrl );
        if( zeroVector )
            nullity = k+1;
        else
            nullity = Min(nullity,k);

        const Real rho_km1_km1 = QR.GetRealPart(k-1,k-1);
        const F rho_km1_k = QR.Get(k-1,k);
        const Real rho_k_k = QR.GetRealPart(k,k); 
        if( ctrl.delta*rho_km1_km1*rho_km1_km1 <=
            rho_k_k*rho_k_k + Abs(rho_km1_k)*Abs(rho_km1_k) )
        {
            ++k;
        }
        else
        {
            ++numSwaps;
            if( ctrl.progress )
                Output("Dropping from k=",k," to ",Max(k-1,1));
            ColSwap( B, k-1, k );
            if( formU )
                ColSwap( U, k-1, k );
            if( formUInv )
                RowSwap( UInv, k-1, k );
            if( k == 1 )
            {
                // We must reinitialize since we keep k=1
                lll::ExpandQR( 0, B, QR, t, d, ctrl.time );
                lll::HouseholderStep( 0, QR, t, d, ctrl.time );
                {
                    auto b0 = B(ALL,IR(0));
                    if( FrobeniusNorm(b0) <= ctrl.zeroTol )
                    {
                        auto QR0 = QR(ALL,IR(0));
                        Zero( b0 );
                        Zero( QR0 );
                        nullity = 1;
                    }
                    else
                        nullity = 0;
                }
            }
            else
            {
                k = k-1; 
            }
        }
    }

    if( ctrl.time )
    {
        Output("  Apply Householder time: ",applyHouseTimer.Total());
        Output("  Round time:             ",roundTimer.Total());
    }

    LLLInfo info;
    info.nullity = nullity;
    info.numSwaps = numSwaps; 
    return info;
}

template<typename F>
LLLInfo BlockedAlg
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  bool formU,
  bool formUInv,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("lll::BlockedAlg"))
    typedef Base<F> Real;
    if( ctrl.time )
    {
        applyHouseTimer.Reset();
        roundTimer.Reset();
        formSInvTimer.Reset();
    }

    const Int m = B.Height();
    const Int n = B.Width();
    const Int minDim = Min(m,n);
    Matrix<F> V, SInv, t;
    Matrix<Real> d;
    Zeros( QR, m, n );
    Zeros( V, m, minDim );
    Zeros( SInv, minDim, minDim );
    Zeros( d, minDim, 1 );
    Zeros( t, minDim, 1 );

    // Perform the first step of Householder reduction
    lll::ExpandBlockQR( 0, B, QR, V, SInv, d, ctrl.time );
    lll::BlockHouseholderStep( 0, QR, V, SInv, t, d, ctrl.time );
    Int nullity = 0;
    {
        auto b0 = B(ALL,IR(0));
        if( FrobeniusNorm(b0) <= ctrl.zeroTol )
        {
            auto QR0 = QR(ALL,IR(0));
            Zero( b0 );
            Zero( QR0 );
            nullity = 1;
        }
    }

    Int k=1, numSwaps=0;
    while( k < n )
    {
        bool zeroVector =
          lll::BlockStep
          ( k, B, U, UInv, QR, V, SInv, t, d, formU, formUInv, ctrl );
        if( zeroVector )
            nullity = k+1;
        else
            nullity = Min(nullity,k);

        const Real rho_km1_km1 = QR.GetRealPart(k-1,k-1);
        const F rho_km1_k = QR.Get(k-1,k);
        const Real rho_k_k = QR.GetRealPart(k,k); 
        if( ctrl.delta*rho_km1_km1*rho_km1_km1 <=
            rho_k_k*rho_k_k + Abs(rho_km1_k)*Abs(rho_km1_k) )
        {
            ++k;
        }
        else
        {
            ++numSwaps;
            if( ctrl.progress )
                Output("Dropping from k=",k," to ",Max(k-1,1));
            ColSwap( B, k-1, k );
            if( formU )
                ColSwap( U, k-1, k );
            if( formUInv )
                RowSwap( UInv, k-1, k );
            if( k == 1 )
            {
                // We must reinitialize since we keep k=1
                lll::ExpandBlockQR( 0, B, QR, V, SInv, d, ctrl.time );
                lll::BlockHouseholderStep( 0, QR, V, SInv, t, d, ctrl.time );
                {
                    auto b0 = B(ALL,IR(0));
                    if( FrobeniusNorm(b0) <= ctrl.zeroTol )
                    {
                        auto QR0 = QR(ALL,IR(0));
                        Zero( b0 );
                        Zero( QR0 );
                        nullity = 1;
                    }
                    else
                        nullity = 0;
                }
            }
            else
            {
                k = k-1; 
            }
        }
    }

    if( ctrl.time )
    {
        Output("  Apply Householder time: ",applyHouseTimer.Total());
        Output("  Form SInv time:         ",formSInvTimer.Total());
        Output("  Round time:             ",roundTimer.Total());
    }

    LLLInfo info;
    info.nullity = nullity;
    info.numSwaps = numSwaps;
    return info;
}

} // namespace lll

template<typename F>
LLLInfo LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    typedef Base<F> Real;
    if( ctrl.delta > Real(1) )
        LogicError("delta is assumed to be at most 1");

    const Int n = B.Width();
    Identity( U, n, n ); 
    Identity( UInv, n, n );

    if( ctrl.presort )
    {
        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> t;
        Matrix<Real> d;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, t, d, Omega, qrCtrl );
        Omega.PermuteCols( B );
        Omega.PermuteCols( U );
        Omega.PermuteRows( UInv );
    }

    const bool useBlocked = false;
    const bool formU = true;
    const bool formUInv = true;
    if( useBlocked )
        return lll::BlockedAlg( B, U, UInv, QR, formU, formUInv, ctrl );
    else
        return lll::UnblockedAlg( B, U, UInv, QR, formU, formUInv, ctrl );
}

template<typename F>
LLLInfo LLL
( Matrix<F>& B,
  Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("LLL"))
    typedef Base<F> Real;
    if( ctrl.delta > Real(1) )
        LogicError("delta is assumed to be at most 1");

    if( ctrl.presort )
    {
        QRCtrl<Real> qrCtrl;
        qrCtrl.smallestFirst = ctrl.smallestFirst;

        auto BCopy = B;
        Matrix<F> t;
        Matrix<Real> d;
        Permutation Omega;
        // TODO: Add support for qr::ProxyHouseholder as well
        El::QR( BCopy, t, d, Omega, qrCtrl );
        Omega.PermuteCols( B );
    }

    const bool useBlocked = false;
    const bool formU = false;
    const bool formUInv = false;
    Matrix<F> U, UInv;
    if( useBlocked )
        return lll::BlockedAlg( B, U, UInv, QR, formU, formUInv, ctrl );
    else
        return lll::UnblockedAlg( B, U, UInv, QR, formU, formUInv, ctrl );
}

template<typename F>
Base<F> LLLDelta
( const Matrix<F>& QR,
  const LLLCtrl<Base<F>>& ctrl )
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
    // for 0 <= k < Min(m,n)-1.
    //
    // TODO: Decide if m < n requires checking the k=m-1 case.
    //
    if( n <= 1 )
        return 1; // the best-possible delta
    Matrix<F> z;
    Real delta = limits::Max<Real>();
    for( Int i=0; i<minDim-1; ++i )
    {
        const Real rho_i_i = R.GetRealPart(i,i);
        if( Abs(rho_i_i) <= ctrl.zeroTol )
            continue;
        const Real rho_i_ip1 = Abs(R.Get(i,i+1));
        const Real rho_ip1_ip1 = R.GetRealPart(i+1,i+1);

        const Real deltaBound =
          (rho_ip1_ip1*rho_ip1_ip1+rho_i_ip1*rho_i_ip1)/(rho_i_i*rho_i_i);
 
        delta = Min(delta,deltaBound);
    }

    // Ensure that
    //
    //    | R(l,k) | <= 0.5 | R(l,l) | for all 0 <= l < Min(m,n) and l < k < n
    //
    // unless a weak reduction was requested in which case k=l+1.
    //
    // TODO: Decide if m < n requires checking the k=m-1 case.
    //
    // NOTE: This does not seem to hold for complex LLL reductions,
    //       so sqrt(2)/2 is used instead.
    const Real bound = 
      ( IsComplex<F>::value ? Real(1)/Sqrt(Real(2)) : Real(1)/Real(2) );
    Real maxRatio = 0;
    if( ctrl.weak )
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R.Get(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                maxRatio = Max(maxRatio,Abs(R.Get(i,i+1)/rho_ii));
        }
    }
    else
    {
        for( Int i=0; i<minDim-1; ++i )
        {
            const F rho_ii = R.Get(i,i);
            if( Abs(rho_ii) <= ctrl.zeroTol )
                continue;
            else
                for( Int j=i+1; j<n; ++j )
                    maxRatio = Max(maxRatio,Abs(R.Get(i,j)/rho_ii));
        }
    }
    if( maxRatio > bound+Pow(limits::Epsilon<Real>(),Real(3)/Real(4)) )
        return 0; // the worst-possible delta

    return delta;
}

#define PROTO(F) \
  template LLLInfo LLL \
  ( Matrix<F>& B, \
    Matrix<F>& QR, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template LLLInfo LLL \
  ( Matrix<F>& B, \
    Matrix<F>& U, \
    Matrix<F>& UInv, \
    Matrix<F>& QR, \
    const LLLCtrl<Base<F>>& ctrl ); \
  template Base<F> LLLDelta \
  ( const Matrix<F>& QR, \
    const LLLCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
