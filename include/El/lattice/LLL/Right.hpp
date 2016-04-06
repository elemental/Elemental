/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_LLL_UNBLOCKED_RIGHT_HPP
#define EL_LATTICE_LLL_UNBLOCKED_RIGHT_HPP

namespace El {
namespace lll {

template<typename F>
void RightGivensStep
( Int k,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d, 
  Matrix<F>& GivensBlock,
  Int GivensFirstCol,
  Int GivensLastCol, 
  bool time )
{        
    if( time )
        formGivensTimer.Start();

    DEBUG_ONLY(CSE cse("lll::RightGivensStep"))
    typedef Base<F> Real;
    const Int m = QR.Height();
    const Int n = QR.Width();
    const Int gBlockSize = GivensBlock.Height();
    const Int gActSize = GivensLastCol - GivensFirstCol + 1;
    const Int colIx = gBlockSize - (GivensLastCol - k) - 1;
    
    if( k == Min(m,n)-1 )
    {
//        if( time )
//            formGivensTimer.Start();
        
        const Real rho_k = QR.GetRealPart(k,k);
        if (rho_k < Real(0))
        {
            QR.Set(k, k, -QR.Get(k,k));
            
            for (Int i=gBlockSize-gActSize-1; i < gBlockSize; i++)
            {
                GivensBlock.Set(i, colIx, -GivensBlock.Get(i,colIx));
            }
        }
        
        if( time )
            formGivensTimer.Stop();
        
        return;
    }
    else if( k >= Min(m,n) )
    {
        if( time )
            formGivensTimer.Stop();
        return;
    }
        
    // Perform one Givens rotation to rows j, j+1
    Real c; F s;
    // G * [x w] = [r u], r, v > 0
    //     [y z]   [0 v]
    F x = QR.Get(k,k); F y = QR.Get(k+1,k);
    F w = QR.Get(k,k+1); F z = QR.Get(k+1,k+1);
    
    lapack::Givens( x, y, &c, &s );
    Matrix<F> G1(2,2);
    c = Sgn(RealPart(x))*Abs(c);
    s = Sgn(RealPart(y))*Abs(s);
    Real sgn = Sgn(RealPart(-Conj(s)*w+c*z));
    
    G1.Set(0,0,c);
    G1.Set(0,1,s);
    G1.Set(1,0,-sgn*Conj(s));
    G1.Set(1,1,sgn*c);
    
    auto RR = QR( IR(k,k+2), IR(k, GivensLastCol+1) );
    Transform2x2Rows(G1, RR, 0, 1);
//    RotateRows( c, s, RR, 0, 1);
    
    auto G = GivensBlock( IR(gBlockSize-gActSize, END), IR(colIx, colIx+2) );
    Transform2x2Cols(G1, G, 0, 1);
//    RotateCols( c, s, G, 0, 1);

    if( time )
        formGivensTimer.Stop();
}

template<typename F>
void RightNegateRow
( Int k,
  Matrix<F>& QR,
  Matrix<F>& GivensBlock,
  Int GivensFirstCol,
  Int GivensLastCol, 
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::RightNegateRow"))
    typedef Base<F> Real;
    const Int n = QR.Width();

    const Real rho_k = QR.GetRealPart(k,k); 
    if (rho_k < Real(0))
    {        
        for (Int i=k; i < n; i++)
        {
            QR.Set(k, i, -QR.Get(k, i));
        }
        if (GivensLastCol > 0)
        {
            const Int gBlockSize = GivensBlock.Height();
            const Int gActSize = GivensLastCol - GivensFirstCol + 1;
            const Int colIx = gBlockSize - (k-GivensFirstCol)-1;
            for (Int i=gBlockSize-gActSize; i < gBlockSize; i++)
            {
                GivensBlock.Set(i, colIx, -GivensBlock.Get(i,colIx));
            }
        }
    }
}

// Put the k'th column of B into the k'th column of QR and then rotate
// said column with the first k-1 (scaled) Householder reflectors.

template<typename F>
void RightExpandQR
( Int k,
        Matrix<F>& QR,
  const Matrix<F>& t,
  const Matrix<Base<F>>& d,
  Int numOrthog,
  Int hPanelStart, 
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::RightExpandQR"))
    typedef Base<F> Real;
    const Int m = QR.Height();
    const Int n = QR.Width();
    const Int minDim = Min(m,n);
          F* QRBuf = QR.Buffer();  
    const F* tBuf = t.LockedBuffer();
    const Base<F>* dBuf = d.LockedBuffer();
    const Int QRLDim = QR.LDim();

    if( time )
        applyHouseTimer.Start();
    for( Int orthog=0; orthog<numOrthog; ++orthog )
    {
        for( Int i=hPanelStart; i<Min(k,minDim); ++i )
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
    }
    if( time )
        applyHouseTimer.Stop();
}

template<typename F>
void RightHouseholderStep
( Int k,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool time )
{
    DEBUG_ONLY(CSE cse("lll::RightHouseholderStep"))
    typedef Base<F> Real;
    const Int m = QR.Height();
    const Int n = QR.Width();
    if( k >= Min(m,n) )
        return;

    if( time )
        houseStepTimer.Start();

    // Perform the next step of Householder reduction
    F* QRBuf = QR.Buffer();
    const Int QRLDim = QR.LDim();
    F& rhokk = QRBuf[k+k*QRLDim]; 
    if( time )
        houseViewTimer.Start();
    auto qr21 = QR( IR(k+1,END), IR(k) );
    if( time )
        houseViewTimer.Stop();
    if( time )
        houseReflectTimer.Start();
    F tau = LeftReflector( rhokk, qr21 );
    if( time )
        houseReflectTimer.Stop();
    t.Set( k, 0, tau );
    if( RealPart(rhokk) < Real(0) )
    {
        d.Set( k, 0, -1 );
        rhokk *= -1;
    }
    else
        d.Set( k, 0, +1 );

    if( time )
        houseStepTimer.Stop();
}

// Return true if the new column is a zero vector
template<typename Z, typename F>
bool RightStep
( Int k,
  Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool formU,
  Int hPanelStart,
  Matrix<F> GivensBlock,
  Int GivensFirstCol,
  Int GivensLastCol,
  bool useHouseholder,
  const LLLCtrl<Base<F>>& ctrl,
  bool& updateCol   )
{
    DEBUG_ONLY(CSE cse("lll::RightStep"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    const Real eps = limits::Epsilon<Real>();

    if( ctrl.time )
        stepTimer.Start();

    Z* BBuf = B.Buffer();
    Z* UBuf = U.Buffer();
    F* QRBuf = QR.Buffer();
    const Int BLDim = B.LDim();
    const Int ULDim = U.LDim();
    const Int QRLDim = QR.LDim();

    while( true ) 
    {
        if( useHouseholder )
            lll::RightExpandQR( k, QR, t, d, ctrl.numOrthog, hPanelStart, ctrl.time);
        else
            lll::RightNegateRow( k, QR, GivensBlock, GivensFirstCol, GivensLastCol, ctrl.time );

        const Base<Z> oldNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
//        if( !limits::IsFinite(oldNorm) )
//            RuntimeError("Encountered an unbounded norm; increase precision");
//        if( oldNorm > Real(1)/eps )
//            RuntimeError("Encountered norm greater than 1/eps, where eps=",eps);

        if( oldNorm <= ctrl.zeroTol )
        {
            for( Int i=0; i<m; ++i )
                BBuf[i+k*BLDim] = 0;
            for( Int i=0; i<m; ++i )
                QRBuf[i+k*QRLDim] = 0;
            if( k < Min(m,n) )
            {
                t.Set( k, 0, Real(2) );
                d.Set( k, 0, Real(1) );
            }
            if( ctrl.time )
                stepTimer.Stop();
            return true;
        }

        if( ctrl.time )
            roundTimer.Start();
        if( ctrl.variant == LLL_WEAK )
        {
            const Real rho_km1_km1 = RealPart(QRBuf[(k-1)+(k-1)*QRLDim]);
            if( rho_km1_km1 > ctrl.zeroTol )
            {
                // TODO: Add while loop?
                F chi = QRBuf[(k-1)+k*QRLDim]/rho_km1_km1;
                if( Abs(RealPart(chi)) > ctrl.eta ||
                    Abs(ImagPart(chi)) > ctrl.eta )
                {
                    chi = Round(chi);
                    blas::Axpy
                    ( k, -chi,
                      &QRBuf[(k-1)*QRLDim], 1,
                      &QRBuf[ k   *QRLDim], 1 );

                    blas::Axpy
                    ( m, -Z(chi),
                      &BBuf[(k-1)*BLDim], 1,
                      &BBuf[ k   *BLDim], 1 );

                    if( formU )
                        blas::Axpy
                        ( n, -Z(chi),
                          &UBuf[(k-1)*ULDim], 1,
                          &UBuf[ k   *ULDim], 1 );
                }
            }
        }
        else
        {
            vector<F> xBuf(k);
            // NOTE: Unless LLL is being aggressively executed in low precision,
            //       this loop should only need to be executed once
            const Int maxSizeReductions = 128;
            for( Int reduce=0; reduce<maxSizeReductions; ++reduce )
            {
                Int numNonzero = 0;
                for( Int i=k-1; i>=0; --i )
                {
                    F chi = QRBuf[i+k*QRLDim]/QRBuf[i+i*QRLDim];
                    if( Abs(RealPart(chi)) > ctrl.eta ||
                        Abs(ImagPart(chi)) > ctrl.eta )
                    {
                        chi = Round(chi);
                        blas::Axpy
                        ( i+1, -chi,
                          &QRBuf[i*QRLDim], 1,
                          &QRBuf[k*QRLDim], 1 );
                        ++numNonzero;
                    }
                    else
                        chi = 0;
                    xBuf[i] = chi;
                }
                
                if( numNonzero == 0 )
                {
                    break;
                }

                updateCol = true;
                
                const float nonzeroRatio = float(numNonzero)/float(k); 
                if( nonzeroRatio >= ctrl.blockingThresh )
                {
                    vector<Z> xzBuf(k);
                    // Need array of type Z
                    for( Int i=0; i<k; ++i)
                        xzBuf[i] = Z(xBuf[i]);
                
                    blas::Gemv
                    ( 'N', m, k,
                      Z(-1), &BBuf[0*BLDim], BLDim,
                             &xzBuf[0],       1,
                      Z(+1), &BBuf[k*BLDim], 1 );
                    if( formU )
                        blas::Gemv
                        ( 'N', n, k,
                          Z(-1), &UBuf[0*ULDim], ULDim,
                                 &xzBuf[0],       1,
                          Z(+1), &UBuf[k*ULDim], 1 );
                }
                else
                {
                    for( Int i=k-1; i>=0; --i )
                    {
                        const Z chi = Z(xBuf[i]);
                        if( chi == Z(0) )
                            continue;
                        blas::Axpy
                        ( m, -chi,
                          &BBuf[i*BLDim], 1,
                          &BBuf[k*BLDim], 1 );
                        if( formU )
                            blas::Axpy
                            ( n, -chi,
                              &UBuf[i*ULDim], 1,
                              &UBuf[k*ULDim], 1 );
                    }
                }
            }
        }
        const Base<Z> newNorm = blas::Nrm2( m, &BBuf[k*BLDim], 1 );
        if( ctrl.time )
            roundTimer.Stop();
//        if( !limits::IsFinite(newNorm) )
//            RuntimeError("Encountered an unbounded norm; increase precision");
//        if( newNorm > Real(1)/eps )
//            RuntimeError("Encountered norm greater than 1/eps, where eps=",eps);

        if( newNorm > ctrl.reorthogTol*oldNorm )
        {
            break;
        }
        else if( ctrl.progress )
            Output
            ("  Reorthogonalizing with k=",k,
             " since oldNorm=",oldNorm," and newNorm=",newNorm);
    }
    if( useHouseholder )
        lll::RightHouseholderStep( k, QR, t, d, ctrl.time );
    if( ctrl.time )
        stepTimer.Stop();
    return false;
}

// NOTE:
// Blocking Givens will screw up linearly dependent columns
// May have to dump in the event of zero columns...keep this in mind

// Consider explicitly returning both Q and R rather than just R (in 'QR')
template<typename Z, typename F>
LLLInfo<Base<F>> RightAlg
( Matrix<Z>& B,
  Matrix<Z>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  bool formU,
  const LLLCtrl<Base<F>>& ctrl)
{
    DEBUG_ONLY(CSE cse("lll::RightAlg"))
    typedef Base<F> Real;
    if( ctrl.time )
    {
        stepTimer.Reset();
        houseStepTimer.Reset();
        houseViewTimer.Reset();
        houseReflectTimer.Reset();
        applyHouseTimer.Reset();
        roundTimer.Reset();
        formQRTimer.Reset();
        applyGivensTimer.Reset();
        copyGivensTimer.Reset();
        formGivensTimer.Reset();
        colNormTimer.Reset();
        LLLTimer.Reset();
        
        LLLTimer.Start();
    }

    // TODO: Make tunable
    const int H_BLOCK_SIZE = 32;
    const Base<F> thresh = Pow(limits::Epsilon<Real>(),Real(0.5));
    if (ctrl.progress)
        cout << "thresh = " << thresh << endl;
    //---------------------------------------
    
    const Int m = B.Height();
    const Int n = B.Width();
    const Int minDim = Min(m,n);
    
    const int G_BLOCK_SIZE = Min(ctrl.givensBlockSize, minDim);
    
    Copy(B, QR);
    Matrix<Base<F>> colNorms;
    Zeros( colNorms, n, 1 );
	
    if( ctrl.time )
        colNormTimer.Start();
        
    // Obtain norms of all columns
    for (Int i=0; i<n; i++)
    {
        auto col = B( ALL, IR(i) );
        colNorms.Set( i, 0, Real(El::FrobeniusNorm(col)) );
    }

    if( ctrl.time )
        colNormTimer.Stop();
    
    Int numSwaps=0;
    Int nullity = 0;
    Int numColFactored = 0;
    if( ctrl.jumpstart && ctrl.startCol > 0 )
    {
        if( QR.Height() != m || QR.Width() != n )
            LogicError
            ("QR was ",QR.Height()," x ",QR.Width()," and should have been ",
             m," x ",n);
        if( t.Height() != minDim || t.Width() != 1 )
            LogicError
            ("t was ",t.Height(),", x ",t.Width()," and should have been ",
             "Min(m,n)=Min(",m,",",n,")=",Min(m,n)," x 1");
        if( d.Height() != minDim || d.Width() != 1 )
            LogicError
            ("d was ",d.Height(),", x ",d.Width()," and should have been ",
             "Min(m,n)=Min(",m,",",n,")=",Min(m,n)," x 1");
        
        auto QR1 = QR( ALL, IR(0, ctrl.startCol));
        auto t1 = t( IR(0, ctrl.startCol), ALL );
        auto d1 = d( IR(0, ctrl.startCol), ALL);
        El::QR(QR1, t1, d1);
        auto QR2 = QR( ALL, IR(ctrl.startCol, END) );
        qr::ApplyQ( LEFT, ADJOINT, QR1, t1, d1, QR2);
        
        numColFactored = ctrl.startCol;     
    }
    else
    {
        Zeros( t, minDim, 1 );
        Zeros( d, minDim, 1 );
        while( true )
        {
            // Perform the first step of Householder reduction
            lll::RightExpandQR( 0, QR, t, d, ctrl.numOrthog, 0, ctrl.time );
            lll::RightHouseholderStep( 0, QR, t, d, ctrl.time );
            if( QR.GetRealPart(0,0) <= ctrl.zeroTol )
            {
                auto b0 = B(ALL,IR(0));
                auto QR0 = QR(ALL,IR(0));
                Zero( b0 );
                Zero( QR0 );
                t.Set( 0, 0, Real(2) );
                d.Set( 0, 0, Real(1) );

                ColSwap( B, 0, (n-1)-nullity );
                ColSwap( QR, 0, (n-1)-nullity );
                if( formU )
                    ColSwap( U, 0, (n-1)-nullity );

                ++nullity;
                ++numSwaps;
            }
            else
            {
                numColFactored = 1;
                break;
            }
            if( nullity >= n )
                break;
        }
    }

    Matrix<F> GivensBlock;
    Identity(GivensBlock, G_BLOCK_SIZE, G_BLOCK_SIZE);
    Int GivensLastCol = -1;
    Int GivensFirstCol = -1;
    
    Int k = ( ctrl.jumpstart ? Max(ctrl.startCol,1) : 1 );
    Int hPanelStart = ( ctrl.jumpstart ? Max(ctrl.startCol,0) : 0 );
    Int hPanelEnd = k;
    bool useHouseholder = true;
    bool updateCol = false;
    while( k < n-nullity )
    {
        if (hPanelEnd - hPanelStart >= H_BLOCK_SIZE)
        {
            // Apply panel of Householder matrices to right
            
            auto H = QR( IR(hPanelStart, END), IR(hPanelStart, hPanelEnd) );
            auto QR2 = QR( IR(hPanelStart, END), IR(k, END) );
            auto t2 = t( IR(hPanelStart, hPanelEnd), ALL );
            auto d2 = d( IR(hPanelStart, hPanelEnd), ALL);
            qr::ApplyQ( LEFT, ADJOINT, H, t2, d2, QR2);
        
            hPanelStart = hPanelEnd;
        }
    
        updateCol = false;
        bool zeroVector = lll::RightStep( k, B, U, QR, t, d, formU, 
            hPanelStart, GivensBlock, GivensFirstCol, GivensLastCol, useHouseholder, ctrl, updateCol );
        
        if( zeroVector )
        {
            ColSwap( B, k, (n-1)-nullity );
            ColSwap( QR, k, (n-1)-nullity );
            if( formU )
                ColSwap( U, k, (n-1)-nullity );
            ++nullity;
            ++numSwaps;
            
            // Update column norms
            colNorms.Set( k, 0, colNorms.Get((n-1)-nullity, 0) );
            colNorms.Set( (n-1)-nullity, 0, 0 );
            
            continue;
        }
        
        if (updateCol)
        {
            if( ctrl.time )
                colNormTimer.Start();
                
            // Update column norm
            auto col = B( ALL, IR(k) );
            colNorms.Set( k, 0, Real(El::FrobeniusNorm(col)) );
            
            auto rCol = QR( IR(0,k+1), IR(k) );
            Base<F> rnorm = El::FrobeniusNorm(rCol);
            
            if( ctrl.time )
                colNormTimer.Stop();
                
            if (ctrl.progress) {
                cout << "k=" << k << ": ||b_k|| = " << colNorms.Get(k, 0) << endl;
                cout << "k=" << k << ": ||r_k|| = " << rnorm << endl;
                cout << "k=" << k << ": (||b_k|| - ||r_k||)/||b_k|| = " << (colNorms.Get(k, 0) - rnorm)/colNorms.Get(k,0) << endl;
            }
            
            if ((colNorms.Get(k, 0) - rnorm)/colNorms.Get(k,0) >= thresh)
            {
                // Dump Givens to the right
                if (GivensLastCol > 0)
                {
                    if( ctrl.time )
                        applyGivensTimer.Start();
                    
                    const Int startIx = G_BLOCK_SIZE-(GivensLastCol - GivensFirstCol) - 1;
                    auto G = GivensBlock(IR(startIx, END), IR(startIx, END));
                    auto QR1 = QR(IR(GivensFirstCol, GivensLastCol+1), IR(GivensLastCol+1, END));
                    
                    if( ctrl.time )
                        copyGivensTimer.Start();
                    
                    Matrix<F> R;
                    Copy(QR1, R);
                    
                    if( ctrl.time )
                        copyGivensTimer.Stop();
                    
                    Gemm(ADJOINT, NORMAL, F(1), G, R, QR1);
                    
                    Identity(GivensBlock, G_BLOCK_SIZE, G_BLOCK_SIZE); // Todo: Change only relevant submatrix to identity.
                    
                    GivensLastCol = k;
                    GivensFirstCol = k-1;
                    
                    if( ctrl.time )
                        applyGivensTimer.Stop();
                }
            
                if (ctrl.progress)
                    cout << "Redoing QR after " << numSwaps  << " swaps" << endl;
                
                if ( ctrl.time )
                    formQRTimer.Start();
                    
                Copy(B, QR);
                El::QR(QR, t, d);
                
                if ( ctrl.time )
                    formQRTimer.Stop();
            }
        }
        
        numColFactored = Max(numColFactored, k+1);
        if (useHouseholder)
            hPanelEnd++;
        
        const Real rho_km1_km1 = QR.GetRealPart(k-1,k-1);
        const F rho_km1_k = QR.Get(k-1,k);
        const Real rho_k_k = ( k >= m ? Real(0) : QR.GetRealPart(k,k) ); 
        
        const Real leftTerm = Sqrt(ctrl.delta)*rho_km1_km1;
        const Real rightTerm = lapack::SafeNorm(rho_k_k,rho_km1_k);
        // NOTE: It is possible that, if delta < 1/2, that rho_k_k could be
        //       zero and the usual Lovasz condition would be satisifed.
        //       For this reason, we explicitly force a pivot if R(k,k) is
        //       deemed to be numerically zero.
        if( leftTerm <= rightTerm && rho_k_k > ctrl.zeroTol )
        {
            ++k;
            if (k >= numColFactored)
            {
                useHouseholder = true;
            }
            
            if (k > GivensLastCol && GivensLastCol > 0)
            {
                if( ctrl.time )
                    applyGivensTimer.Start();
                    
                const Int startIx = G_BLOCK_SIZE-(GivensLastCol - GivensFirstCol) - 1;
                auto G = GivensBlock(IR(startIx, END), IR(startIx, END));
                auto QR1 = QR(IR(GivensFirstCol, GivensLastCol+1), IR(GivensLastCol+1, END));
                if( ctrl.time )
                    copyGivensTimer.Start();
                Matrix<F> R;
                Copy(QR1, R);
                if( ctrl.time )
                    copyGivensTimer.Stop();
                
                Gemm(ADJOINT, NORMAL, F(1), G, R, QR1);
                
                Identity(GivensBlock, G_BLOCK_SIZE, G_BLOCK_SIZE); // Todo: Change only relevant submatrix to identity.
                GivensLastCol = -1;
                GivensFirstCol = -1;
                
                if( ctrl.time )
                    applyGivensTimer.Stop();
            }
        }
        else
        {
            if (useHouseholder)
            {
                // Apply panel of Householder matrices to right
                auto H = QR( IR(hPanelStart, END), IR(hPanelStart, hPanelEnd) );
                auto QR2 = QR( IR(hPanelStart, END), IR(k+1, END) );
                auto t2 = t( IR(hPanelStart, hPanelEnd), ALL );
                auto d2 = d( IR(hPanelStart, hPanelEnd), ALL);
                qr::ApplyQ( LEFT, ADJOINT, H, t2, d2, QR2);
                        
                hPanelStart = hPanelEnd;
                
                // Switch to Givens regime
                useHouseholder = false;
                
                GivensLastCol = k;
                GivensFirstCol = k-1;
            }
        
            if (k-1 <= GivensLastCol - G_BLOCK_SIZE)
            {
                if( ctrl.time )
                    applyGivensTimer.Start();
                    
                // Dump Givens rotations right of GivensLastCol
                const Int startIx = G_BLOCK_SIZE-(GivensLastCol - GivensFirstCol) - 1;
                auto G = GivensBlock(IR(startIx, END), IR(startIx, END));
                auto QR1 = QR(IR(GivensFirstCol, GivensLastCol+1), IR(GivensLastCol+1, END));
                if( ctrl.time )
                    copyGivensTimer.Start();
                Matrix<F> R;
                Copy(QR1, R);
                if( ctrl.time )
                    copyGivensTimer.Stop();
                
                Gemm(ADJOINT, NORMAL, F(1), G, R, QR1);
                
                Identity(GivensBlock, G_BLOCK_SIZE, G_BLOCK_SIZE); // Todo: Change only relevant submatrix to identity.
                
                GivensLastCol = k;
                GivensFirstCol = k-1;
                
                if( ctrl.time )
                    applyGivensTimer.Stop();
            }
        
            if (GivensLastCol < 0)
            {
                GivensLastCol = k;
                GivensFirstCol = k-1;
            }
        
            ++numSwaps;
            if( ctrl.progress )
            {
                if( rho_k_k <= ctrl.zeroTol )
                    Output("Dropping from k=",k," because R(k,k) ~= 0");
                else
                    Output
                    ("Dropping from k=",k," to ",Max(k-1,1),
                     " since sqrt(delta)*R(k-1,k-1)=",leftTerm," > ",rightTerm);
            }

            
            ColSwap( B, k-1, k );
            ColSwap( QR, k-1, k);
            QR.Set(k,k,F(0)); // This entry doesn't exist in R
            
            if( ctrl.time )
                colNormTimer.Start();
            
            // Update column norms
            Base<F> tmp = colNorms.Get(k, 0);
            colNorms.Set( k, 0, colNorms.Get(k-1, 0) );
            colNorms.Set( k-1, 0, tmp );
            
            if( ctrl.time )
                colNormTimer.Stop();
            
            if( formU )
                ColSwap( U, k-1, k );
                
            if( k == 1 )
            {
                GivensFirstCol = 0;
            
                while( true )
                {
                    // We must reinitialize since we keep k=1
                    lll::RightGivensStep( 0, QR, t, d, GivensBlock, GivensFirstCol, GivensLastCol, ctrl.time );
                    
                    if( QR.GetRealPart(0,0) <= ctrl.zeroTol )
                    {
                        auto b0 = B(ALL,IR(0));
                        auto QR0 = QR(ALL,IR(0));
                        Zero( b0 );
                        Zero( QR0 );
                        t.Set( 0, 0, Real(2) );
                        d.Set( 0, 0, Real(1) );

                        ColSwap( B, 0, (n-1)-nullity );
                        ColSwap( QR, 0, (n-1)-nullity );
                        if( formU )
                            ColSwap( U, 0, (n-1)-nullity );
                       
                        ++nullity;
                        ++numSwaps;
                    }
                    else
                        break;
                    if( nullity >= n )
                        break;
                }
            }
            else
            {
                k = k-1;
                // Update the perceived size of the GivensBlock matrix
                GivensFirstCol = Min(k, GivensFirstCol);
                // Immediately correct the swap in the R factor
                lll::RightGivensStep( k, QR, t, d, GivensBlock, GivensFirstCol, GivensLastCol, ctrl.time );
            }
        }
    }

    // TODO: Need to apply H to Q before exiting (if no swaps at end)?
    
    if( ctrl.time )
    {
        LLLTimer.Stop();
    
        Output("  Step time:              ",stepTimer.Total());
        Output("    Householder step time:  ",houseStepTimer.Total());
        Output("      view time:              ",houseViewTimer.Total());
        Output("      reflect time:           ",houseReflectTimer.Total());
        Output("    Apply Householder time: ",applyHouseTimer.Total());
        Output("    Round time:             ",roundTimer.Total());
        Output("  Form QR time:           ",formQRTimer.Total());
        Output("  Apply Givens time:      ",applyGivensTimer.Total());  
        Output("  Copy Givens time:       ",copyGivensTimer.Total());
        Output("  Form Givens time:       ",formGivensTimer.Total());
        Output("  ColNorm time:           ",colNormTimer.Total());
        Output("Total LLL time:           ",LLLTimer.Total());
    }

    std::pair<Real,Real> achieved = lll::Achieved(QR,ctrl);
    Real logVol = lll::LogVolume(QR);

    LLLInfo<Base<F>> info;
    info.delta = achieved.first;
    info.eta = achieved.second;
    info.rank = n-nullity;
    info.nullity = nullity;
    info.numSwaps = numSwaps;
    info.logVol = logVol;

    return info;
}

} // namespace lll
} // namespace El

#endif // ifndef EL_LATTICE_LLL_UNBLOCKED_HPP
