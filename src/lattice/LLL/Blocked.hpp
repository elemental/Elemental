/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_LLL_BLOCKED_HPP
#define EL_LATTICE_LLL_BLOCKED_HPP

namespace El {
namespace lll {

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
    const Real eps = limits::Epsilon<Real>();

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
        if( !limits::IsFinite(oldNorm) )
            RuntimeError("Encountered an unbounded norm; increase precision");
        if( oldNorm > Real(1)/eps )
            RuntimeError("Encountered norm greater than 1/eps, where eps=",eps);
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
                F chi = QRBuf[(k-1)+k*QRLDim]/QRBuf[(k-1)+(k-1)*QRLDim];
                if( Abs(RealPart(chi)) > ctrl.eta ||
                    Abs(ImagPart(chi)) > ctrl.eta )
                {
                    chi = Round(chi);
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
                else
                    chi = 0;
            }
        }
        else
        {
            vector<F> xBuf(k);
            for( Int i=k-1; i>=0; --i )
            {
                if( Abs(QRBuf[i+i*QRLDim]) <= ctrl.zeroTol )
                {
                    xBuf[i] = 0;
                    continue;
                }
                F chi = QRBuf[i+k*QRLDim]/QRBuf[i+i*QRLDim];
                if( Abs(RealPart(chi)) > ctrl.eta ||
                    Abs(ImagPart(chi)) > ctrl.eta )
                {
                    chi = Round(chi);
                    blas::Axpy
                    ( i+1, -chi,
                      &QRBuf[i*QRLDim], 1,
                      &QRBuf[k*QRLDim], 1 );
                }
                else
                    chi = 0;
                xBuf[i] = chi;
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
        if( !limits::IsFinite(newNorm) )
            RuntimeError("Encountered an unbounded norm; increase precision");
        if( newNorm > Real(1)/eps )
            RuntimeError("Encountered norm greater than 1/eps, where eps=",eps);

        if( newNorm > ctrl.reorthogTol*oldNorm )
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

// Consider explicitly returning both Q and R rather than just R (in 'QR')
template<typename F>
LLLInfo<Base<F>> BlockedAlg
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

        const Real leftTerm = Sqrt(ctrl.delta)*rho_km1_km1;
        const Real rightTerm = lapack::SafeNorm(rho_k_k,rho_km1_k);
        if( leftTerm <= rightTerm )
        {
            ++k;
        }
        else
        {
            ++numSwaps;
            if( ctrl.progress )
                Output("Dropping from k=",k," to ",Max(k-1,1)," since sqrt(delta)*R(k-1,k-1)=",leftTerm," > ",rightTerm);
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

    // Force R to be upper-trapezoidal
    MakeTrapezoidal( UPPER, QR );

    std::pair<Real,Real> achieved = lll::Achieved(QR,ctrl);
    Real logAbsDet = lll::LogAbsDet(QR);

    LLLInfo<Base<F>> info;
    info.delta = achieved.first;
    info.eta = achieved.second;
    info.rank = n-nullity;
    info.nullity = nullity;
    info.numSwaps = numSwaps;
    info.logAbsDet = logAbsDet;
    return info;
}

} // namespace lll
} // namespace El

#endif // ifndef EL_LATTICE_LLL_BLOCKED_HPP
