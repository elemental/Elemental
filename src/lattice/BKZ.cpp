/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace bkz {

template<typename F>
bool TrivialCoordinates( const Matrix<F>& v )
{
    DEBUG_ONLY(CSE cse("bkz::TrivialCoordinates"))
    const Int n = v.Height();    
    if( n == 0 )
        LogicError("Invalid coordinate length");
    if( v.Get(0,0) != F(1) )
        return false;
    for( Int i=1; i<n; ++i )
        if( v.Get(i,0) != F(0) )
            return false;
    return true;
}

} // namespace bkz

template<typename F>
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZWithQ"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    if( ctrl.blocksize < 2 )
        LogicError("BKZ requires a blocksize of at least 2");
    if( m < n )
        LogicError("If height(B) < width(B), call LLL first");

    const bool progress = true;

    Int numSwaps=0;
    auto lllInfo = LLLWithQ( B, U, UInv, QR, t, d, ctrl.lllCtrl );
    numSwaps = lllInfo.numSwaps;

    Int z=0, j=-1; 
    Matrix<F> BTmp, UTmp, UInvTmp, QRTmp, tTmp;
    Matrix<Base<F>> dTmp;
    Int numEnums=0, numEnumFailures=0;
    while( z < n-1 ) 
    {
        j = Mod(j+1,n);
        const Int k = Min(j+ctrl.blocksize-1,n-1);
        const Int h = Min(k+1,n-1); 
        
        Matrix<F> v;
        {
            auto BEnum = B( ALL, IR(j,k+1) );
            auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
            ShortestVectorEnumeration( BEnum, QREnum, v, ctrl.probabalistic );
            ++numEnums;
        }

        if( bkz::TrivialCoordinates(v) )        
        {
            if( progress )
                Output("  Trivial enumeration for j=",j,", z=",z);
            ++z;
            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = h-1;
            const auto subInd = IR(0,h+1);
            auto BSub = B( ALL, subInd );
            auto USub = U( subInd, subInd );
            auto UInvSub = UInv( subInd, subInd ); 
            auto QRSub = QR( ALL, subInd );
            auto tSub = t( subInd, ALL );
            auto dSub = d( subInd, ALL );
            lllInfo =
              LLLWithQ( BSub, USub, UInvSub, QRSub, tSub, dSub, subLLLCtrl );
            numSwaps += lllInfo.numSwaps;
        }
        else
        {
            if( progress )
                Output("  Nontrivial enumeration for j=",j,", z=",z);
            ++numEnumFailures;
            z = 0;
            Matrix<F> bNew;
            Zeros( bNew, m, 1 );
            Gemv( NORMAL, F(1), B(ALL,IR(j,k+1)), v, bNew );

            // The following code looks rather complicated but is due to 
            // running LLL on the extended basis
            //
            //   (b_0,b_1,...,b_{j-1},\sum_{i=j}^{k} v_i,b_j,...,b_h)
            //
            // using the fact that, if
            //
            //   [B_L,b,B_R] [a,[G_T;g_M;G_R]] = [0,\tilde{B}],
            // 
            // then
            //
            //   [B_L,B_R] U = \tilde{B},
            //
            // where U = [G_T;G_R] + b g_M. Furthermore, if
            //
            //   [B_L,b,B_R] = [0,\tilde{B}] [q; U^{-1}_L, r, U^{-1}_R],
            //
            // then
            //
            //   [B_L,B_R] = \tilde{B} [U^{-1}_L, U^{-1}_R].
            //
            // The initializations of the expanded U and UInv are simpler to
            // derive.
 
            BTmp.Resize( m, h+2 );
            {
                auto BL = B( ALL, IR(0,j) );
                auto BR = B( ALL, IR(j,h+1) );
                auto BTmpL = BTmp( ALL, IR(0,j)     );
                auto bTmpM = BTmp( ALL, IR(j)       );
                auto BTmpR = BTmp( ALL, IR(j+1,h+2) );
                BTmpL = BL;
                bTmpM = bNew;
                BTmpR = BR;
            }
            Zeros( UTmp, h+2, h+2 );
            Zeros( UInvTmp, h+2, h+2 );
            {
                auto UTL = U( IR(0,j), IR(0,j) );
                auto UTR = U( IR(0,j), IR(j,h+1) );
                auto UBL = U( IR(j,h+1), IR(0,j) );
                auto UBR = U( IR(j,h+1), IR(j,h+1) );
                auto UTmpTL = UTmp( IR(0,j), IR(0,j) );
                auto UTmpTR = UTmp( IR(0,j), IR(j+1,h+2) );
                auto UTmpBL = UTmp( IR(j+1,h+2), IR(0,j) );
                auto UTmpBR = UTmp( IR(j+1,h+2), IR(j+1,h+2) );
                UTmpTL = UTL;
                UTmpTR = UTR;
                UTmpBL = UBL;
                UTmpBR = UBR;
                UTmp.Set( j, j, F(1) );

                auto UInvTL = UInv( IR(0,j), IR(0,j) );
                auto UInvTR = UInv( IR(0,j), IR(j,h+1) );
                auto UInvBL = UInv( IR(j,h+1), IR(0,j) );
                auto UInvBR = UInv( IR(j,h+1), IR(j,h+1) );
                auto UInvTmpTL = UInvTmp( IR(0,j), IR(0,j) );
                auto UInvTmpTR = UInvTmp( IR(0,j), IR(j+1,h+2) );
                auto UInvTmpBL = UInvTmp( IR(j+1,h+2), IR(0,j) );
                auto UInvTmpBR = UInvTmp( IR(j+1,h+2), IR(j+1,h+2) );
                UInvTmpTL = UInvTL;
                UInvTmpTR = UInvTR;
                UInvTmpBL = UInvBL;
                UInvTmpBR = UInvBR;
                UInvTmp.Set( j, j, F(1) );
            }
            QRTmp.Resize( m, h+2 );
            {
                auto QRL = QR( ALL, IR(0,j) );
                auto QRTmpL = QRTmp( ALL, IR(0,j) );
                QRTmpL = QRL;
            }
            tTmp.Resize( Min(m,h+2), 1 );
            dTmp.Resize( Min(m,h+2), 1 );
            {
                auto tT = t( IR(0,j), ALL );
                auto tTmpT = tTmp( IR(0,j), ALL );
                tTmpT = tT;

                auto dT = d( IR(0,j), ALL );
                auto dTmpT = dTmp( IR(0,j), ALL );
                dTmpT = dT;
            }

            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = j;
            lllInfo =
              LLLWithQ( BTmp, UTmp, UInvTmp, QRTmp, tTmp, dTmp, subLLLCtrl );
            numSwaps += lllInfo.numSwaps;
            {
                // The last column of BTmp should be all zeros now
                auto BL = B( ALL, IR(0,h+1) );
                auto BTmpL = BTmp( ALL, IR(0,h+1) );
                BL = BTmpL;
            }
            {
                auto USubT = U( IR(0,j), IR(0,h+1) );
                auto USubB = U( IR(j,h+1), IR(0,h+1) );
                auto UTmpT = UTmp( IR(0,j), IR(0,h+1) );
                auto uTmpM = UTmp( IR(j), IR(0,h+1) );
                auto UTmpB = UTmp( IR(j+1,h+2), IR(0,h+1) );
                USubT = UTmpT;
                USubB = UTmpB;
                auto USub_k = U( IR(j,k+1), IR(0,h+1) );
                Geru( F(1), v, uTmpM, USub_k );

                auto UInvSubL = UInv( IR(0,h+1), IR(0,j) );
                auto UInvSubR = UInv( IR(0,h+1), IR(j,h+1) );
                auto UInvTmpL = UInvTmp( IR(0,h+1), IR(0,j) ); 
                auto uInvTmpM = UInvTmp( IR(0,h+1), IR(j) );
                auto UInvTmpR = UInvTmp( IR(0,h+1), IR(j+1,h+2) );
                UInvSubL = UInvTmpL;
                UInvSubR = UInvTmpR;
            }
            // Returning the QR factorization doesn't work without explicitly
            // forming Q due to the first column of BTmp being removed
            subLLLCtrl.startCol = 0;
            lllInfo = LLLWithQ( B, U, UInv, QR, t, d, ctrl.lllCtrl );
        }

        if( ctrl.earlyAbort && numEnums >= ctrl.numEnumsBeforeAbort )
            break;
    }

    // Perform a final pass to get the full LLL info
    // NOTE: This could be replaced in favor of manually computing the 
    //       returned LLL info but should be cheap relative to the above BKZ
    lllInfo = LLLWithQ( B, U, UInv, QR, t, d, ctrl.lllCtrl );
    numSwaps += lllInfo.numSwaps;

    BKZInfo<Real> info;
    info.delta = lllInfo.delta;
    info.eta = lllInfo.eta;
    info.rank = lllInfo.rank;
    info.nullity = lllInfo.nullity;
    info.numSwaps = numSwaps;
    info.numEnums = numEnums;
    info.numEnumFailures = numEnumFailures;
    info.logVol = lllInfo.logVol;
    return info;
}

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZ"))
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = BKZWithQ( B, U, UInv, R, t, d, ctrl );
    MakeTrapezoidal( UPPER, R );
    return info;
}

template<typename F>
BKZInfo<Base<F>>
BKZWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZWithQ"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    if( ctrl.blocksize < 2 )
        LogicError("BKZ requires a blocksize of at least 2");
    if( m < n )
        LogicError("If height(B) < width(B), call LLL first");

    const bool progress = true;

    Int numSwaps=0;
    auto lllInfo = LLLWithQ( B, QR, t, d, ctrl.lllCtrl );
    numSwaps = lllInfo.numSwaps;

    Int z=0, j=-1; 
    Matrix<F> BTmp, QRTmp, tTmp;
    Matrix<Base<F>> dTmp;
    Int numEnums=0, numEnumFailures=0;
    while( z < n-1 ) 
    {
        j = Mod(j+1,n);
        const Int k = Min(j+ctrl.blocksize-1,n-1);
        const Int h = Min(k+1,n-1); 

        Matrix<F> v;
        {
            auto BEnum = B( ALL, IR(j,k+1) );
            auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
            ShortestVectorEnumeration( BEnum, QREnum, v, ctrl.probabalistic );
            ++numEnums;
        }

        if( bkz::TrivialCoordinates(v) )        
        {
            if( progress )
                Output("  Trivial enumeration for j=",j,", z=",z);
            ++z;
            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = h-1;
            const auto subInd = IR(0,h+1);
            auto BSub = B( ALL, subInd );
            auto QRSub = QR( ALL, subInd );
            auto tSub = t( subInd, ALL );
            auto dSub = d( subInd, ALL );
            lllInfo = LLLWithQ( BSub, QRSub, tSub, dSub, subLLLCtrl );
            numSwaps += lllInfo.numSwaps;
        }
        else
        {
            if( progress )
                Output("  Nontrivial enumeration for j=",j,", z=",z);
            ++numEnumFailures;
            z = 0;
            Matrix<F> bNew;
            Zeros( bNew, m, 1 );
            Gemv( NORMAL, F(1), B(ALL,IR(j,k+1)), v, bNew );

            BTmp.Resize( m, h+2 );
            {
                auto BL = B( ALL, IR(0,j) );
                auto BR = B( ALL, IR(j,h+1) );
                auto BTmpL = BTmp( ALL, IR(0,j)     );
                auto bTmpM = BTmp( ALL, IR(j)       );
                auto BTmpR = BTmp( ALL, IR(j+1,h+2) );
                BTmpL = BL;
                bTmpM = bNew;
                BTmpR = BR;
            }
            QRTmp.Resize( m, h+2 );
            {
                auto QRL = QR( ALL, IR(0,j) );
                auto QRTmpL = QRTmp( ALL, IR(0,j) );
                QRTmpL = QRL;
            }
            tTmp.Resize( Min(m,h+2), 1 );
            dTmp.Resize( Min(m,h+2), 1 );
            {
                auto tT = t( IR(0,j), ALL );
                auto tTmpT = tTmp( IR(0,j), ALL );
                tTmpT = tT;

                auto dT = d( IR(0,j), ALL );
                auto dTmpT = dTmp( IR(0,j), ALL );
                dTmpT = dT;
            }

            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = j;
            auto lllInfo = LLLWithQ( BTmp, QRTmp, tTmp, dTmp, subLLLCtrl );
            numSwaps += lllInfo.numSwaps;
            {
                // The last column of BTmp should be all zeros now
                auto BL = B( ALL, IR(0,h+1) );
                auto BTmpL = BTmp( ALL, IR(0,h+1) );
                BL = BTmpL;
            }
            // Returning the QR factorization doesn't work without explicitly
            // forming Q due to the first column of BTmp being removed
            subLLLCtrl.startCol = 0;
            lllInfo = LLLWithQ( B, QR, t, d, ctrl.lllCtrl );
        }

        if( ctrl.earlyAbort && numEnums >= ctrl.numEnumsBeforeAbort )
            break;
    }

    // Perform a final pass to get the full LLL info
    // NOTE: This could be replaced in favor of manually computing the 
    //       returned LLL info but should be cheap relative to the above BKZ
    lllInfo = LLLWithQ( B, QR, t, d, ctrl.lllCtrl );
    numSwaps += lllInfo.numSwaps;

    BKZInfo<Real> info;
    info.delta = lllInfo.delta;
    info.eta = lllInfo.eta;
    info.rank = lllInfo.rank;
    info.nullity = lllInfo.nullity;
    info.numSwaps = numSwaps;
    info.numEnums = numEnums;
    info.numEnumFailures = numEnumFailures;
    info.logVol = lllInfo.logVol;
    return info;
}

template<typename F>
BKZInfo<Base<F>>
BKZ
( Matrix<F>& B,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZ"))
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = BKZWithQ( B, R, t, d, ctrl );
    MakeTrapezoidal( UPPER, R );
    return info;
}

// Emulate the flavor of quicksort/mergesort by recursively splitting the 
// vectors in half, applying BKZ to each half, and merging the halves
// by running BKZ on the interwoven reduced basis vectors
// (notice that this should allow the highest levels to often run at a lower
//  precision since the reduced basis vectors are likely to be much smaller, 
//  especially with SVP challenge lattices).
//
// C.f. The analogue of Lehmer's version of Euclid's algorithm that Schnorr
// mentions at the end of "Progress on BKZ and Lattice Reduction".
//
// NOTE: Until Complex<BigFloat> exists, we must have different implementations
//       for real and complex F

// TODO: Provide a way to display when changing precision without showing
//       all of the backtracks and deep insertions.

namespace bkz {

template<typename F,typename RealLower>
BKZInfo<RealLower>
LowerPrecisionMerge
( const Matrix<F>& CL,
  const Matrix<F>& CR,
        Matrix<F>& B,
        Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("bkz::LowerPrecisionMerge"))
    typedef ConvertBase<F,RealLower> FLower;
    const string typeString = TypeName<RealLower>();

    const Int n = B.Width();
    const Int firstHalf = n-(n/2);

    if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        Output("  Dropping to " + typeString);
    Matrix<FLower> BLower;
    BLower.Resize( B.Height(), n );
    // Interleave CL and CR to reform B before running BKZ again
    // NOTE: This does not seem to make a substantial difference
    for( Int jSub=0; jSub<n/2; ++jSub )
    {
        auto cl = CL( ALL, IR(jSub) );
        auto cr = CR( ALL, IR(jSub) );
        auto bl = BLower( ALL, IR(2*jSub) );
        auto br = BLower( ALL, IR(2*jSub+1) );
        Copy( cl, bl );
        Copy( cr, br );
    }
    if( firstHalf > n/2 )
    {
        auto cl = CL( ALL, IR(firstHalf-1) );
        auto bl = BLower( ALL, IR(n-1) );
        Copy( cl, bl );
    }

    BKZCtrl<RealLower> ctrlLower( ctrl );
    RealLower eps = limits::Epsilon<RealLower>();
    RealLower minEta = RealLower(1)/RealLower(2)+Pow(eps,RealLower(0.9));
    if( ctrlLower.lllCtrl.eta < minEta )
    {
        Output("eps=",eps);
        Output("ctrlLower.lllCtrl.eta=",ctrlLower.lllCtrl.eta);
        Output("Max(",RealLower(ctrl.lllCtrl.eta),",",minEta,")=",Max(RealLower(ctrl.lllCtrl.eta),minEta));
        ctrlLower.lllCtrl.eta = minEta;
        Output("ctrlLower.lllCtrl.eta new=",ctrlLower.lllCtrl.eta);
    }
    Timer timer;
    Matrix<FLower> RLower;
    if( ctrl.lllCtrl.time )
        timer.Start();
    auto infoLower = BKZ( BLower, RLower, ctrlLower );
    if( ctrl.lllCtrl.time )
        Output("  " + typeString + " BKZ took ",timer.Stop()," seconds");
    Copy( BLower, B );
    Copy( RLower, R );
    return infoLower;
}

template<typename Real>
BKZInfo<Real>
RecursiveHelper
( Matrix<Real>& B,
  Matrix<Real>& R,
  Int numShuffles,
  Int cutoff,
  const BKZCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bkz::RecursiveHelper"))
    typedef Real F;
    const Int n = B.Width();
    if( n < cutoff )
        return BKZ( B, R, ctrl );
    Timer timer;

    // Reduce the entire matrix with LLL before attempting BKZ.
    // Deep reductions should probably not be used due to the expense.
    BKZInfo<Real> info;
    info.numSwaps = 0;
    auto lllInfo = RecursiveLLL( B, R, cutoff, ctrl.lllCtrl );
    info.numSwaps += lllInfo.numSwaps;

    const Real BOneNorm = OneNorm( B );
    Output("|| B ||_1 = ",BOneNorm);

    for( Int shuffle=0; shuffle<=numShuffles; ++shuffle )
    {
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
            Output("Shuffle=",shuffle);
        auto C( B ); 

        const Int firstHalf = n-(n/2);
        auto CL = C( ALL, IR(0,firstHalf) );
        auto CR = C( ALL, IR(firstHalf,n) );

        double leftTime, rightTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        auto leftInfo = RecursiveBKZ( CL, cutoff, ctrl ); 
        if( ctrl.lllCtrl.time )
        {
            leftTime = timer.Stop(); 
            timer.Start();
        }
        auto rightInfo = RecursiveBKZ( CR, cutoff, ctrl );
        if( ctrl.lllCtrl.time )
            rightTime = timer.Stop();
        info.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("n=",n);
            Output("  left swaps=",leftInfo.numSwaps);
            Output("  right swaps=",rightInfo.numSwaps);
        }
        if( ctrl.lllCtrl.time )
        {
            Output("  left time:  ",leftTime," seconds");
            Output("  right time: ",rightTime," seconds");
        }
        const Real CLOneNorm = OneNorm( CL );
        const Real CROneNorm = OneNorm( CR );
        const Real CLMaxNorm = MaxNorm( CL );
        const Real CRMaxNorm = MaxNorm( CR );
        if( ctrl.lllCtrl.progress )
        {
            Output("  || C_L ||_1 = ",CLOneNorm);
            Output("  || C_R ||_1 = ",CROneNorm);
            Output("  || C_L ||_max = ",CLMaxNorm);
            Output("  || C_R ||_max = ",CRMaxNorm);
        }

        const Real COneNorm = Max(CLOneNorm,CROneNorm);
        const Real fudge = 1.5; // TODO: Make tunable
        const Int neededPrec = Int(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( PrecisionIsGreater<Real,float>::value && neededPrec <= 24 )
        {
            try
            {
                info = LowerPrecisionMerge<F,float>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            PrecisionIsGreater<Real,double>::value && neededPrec <= 53 )
        {
            try
            {
                info = LowerPrecisionMerge<F,double>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#ifdef EL_HAVE_QUAD
        if( !succeeded &&
            PrecisionIsGreater<Real,Quad>::value && neededPrec <= 113 )
        {
            try
            {
                info = LowerPrecisionMerge<F,Quad>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#endif
#ifdef EL_HAVE_MPC
        // Only move down to a lower-precision MPFR type if the jump is
        // substantial. The current value has been naively chosen.
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpc::Precision();
        if( !succeeded && neededPrec <= inputPrec-minPrecDiff )
        {
            mpc::SetPrecision( neededPrec );
            try {
                info = LowerPrecisionMerge<F,BigFloat>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
            mpc::SetPrecision( inputPrec );
        }
#endif

        if( !succeeded )
        {
            // Interleave CL and CR to reform B before running BKZ again
            for( Int jSub=0; jSub<n/2; ++jSub )
            {
                auto cl = CL( ALL, IR(jSub) );
                auto cr = CR( ALL, IR(jSub) ); 
                auto bl = B( ALL, IR(2*jSub) );
                auto br = B( ALL, IR(2*jSub+1) );
                bl = cl;
                br = cr;
            }
            if( firstHalf > n/2 )
            {
                auto cl = CL( ALL, IR(firstHalf-1) );
                auto bl = B( ALL, IR(n-1) ); 
                bl = cl;
            }
            
            info = BKZ( B, R, ctrl );
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

// Same as the above, but with the Complex<BigFloat> datatype avoided
template<typename Real>
BKZInfo<Real>
RecursiveHelper
( Matrix<Complex<Real>>& B,
  Matrix<Complex<Real>>& R,
  Int numShuffles,
  Int cutoff,
  const BKZCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bkz::RecursiveHelper"))
    typedef Complex<Real> F;
    const Int n = B.Width();
    if( n < cutoff )
        return BKZ( B, R, ctrl );
    Timer timer;

    // Reduce the entire matrix with LLL before attempting BKZ.
    // Deep reductions should probably not be used due to the expense.
    BKZInfo<Real> info;
    info.numSwaps = 0;
    auto lllInfo = RecursiveLLL( B, R, cutoff, ctrl.lllCtrl );
    info.numSwaps += lllInfo.numSwaps;

    const Real BOneNorm = OneNorm( B );
    Output("|| B ||_1 = ",BOneNorm);

    for( Int shuffle=0; shuffle<=numShuffles; ++shuffle )
    {
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
            Output("Shuffle=",shuffle);
        auto C( B ); 

        const Int firstHalf = n-(n/2);
        auto CL = C( ALL, IR(0,firstHalf) );
        auto CR = C( ALL, IR(firstHalf,n) );

        double leftTime, rightTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        auto leftInfo = RecursiveBKZ( CL, cutoff, ctrl ); 
        if( ctrl.lllCtrl.time )
        {
            leftTime = timer.Stop(); 
            timer.Start();
        }
        auto rightInfo = RecursiveBKZ( CR, cutoff, ctrl );
        if( ctrl.lllCtrl.time )
            rightTime = timer.Stop();
        info.numSwaps += leftInfo.numSwaps + rightInfo.numSwaps;
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("n=",n);
            Output("  left swaps=",leftInfo.numSwaps);
            Output("  right swaps=",rightInfo.numSwaps);
        }
        if( ctrl.lllCtrl.time )
        {
            Output("  left time:  ",leftTime," seconds");
            Output("  right time: ",rightTime," seconds");
        }
        const Real CLOneNorm = OneNorm( CL );
        const Real CROneNorm = OneNorm( CR );
        const Real CLMaxNorm = MaxNorm( CL );
        const Real CRMaxNorm = MaxNorm( CR );
        if( ctrl.lllCtrl.progress )
        {
            Output("  || C_L ||_1 = ",CLOneNorm);
            Output("  || C_R ||_1 = ",CROneNorm);
            Output("  || C_L ||_max = ",CLMaxNorm);
            Output("  || C_R ||_max = ",CRMaxNorm);
        }

        const Real COneNorm = Max(CLOneNorm,CROneNorm);
        const Real fudge = 1.5; // TODO: Make tunable
        const Int neededPrec = Int(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( PrecisionIsGreater<Real,float>::value && neededPrec <= 24 )
        {
            try
            {
                info = LowerPrecisionMerge<F,float>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            PrecisionIsGreater<Real,double>::value && neededPrec <= 53 )
        {
            try
            {
                info = LowerPrecisionMerge<F,double>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
#ifdef EL_HAVE_QUAD
        if( !succeeded &&
            PrecisionIsGreater<Real,Quad>::value && neededPrec <= 113 )
        {
            try
            {
                info = LowerPrecisionMerge<F,Quad>( CL, CR, B, R, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
#endif

        if( !succeeded )
        {
            // Interleave CL and CR to reform B before running BKZ again
            for( Int jSub=0; jSub<n/2; ++jSub )
            {
                auto cl = CL( ALL, IR(jSub) );
                auto cr = CR( ALL, IR(jSub) ); 
                auto bl = B( ALL, IR(2*jSub) );
                auto br = B( ALL, IR(2*jSub+1) );
                bl = cl;
                br = cr;
            }
            if( firstHalf > n/2 )
            {
                auto cl = CL( ALL, IR(firstHalf-1) );
                auto bl = B( ALL, IR(n-1) ); 
                bl = cl;
            }
            
            info = BKZ( B, R, ctrl );
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

} // namespace bkz

template<typename F>
BKZInfo<Base<F>>
RecursiveBKZ
( Matrix<F>& B,
  Int cutoff,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveBKZ"))
    Matrix<F> R;
    return RecursiveBKZ( B, R, cutoff, ctrl );
}

template<typename F>
BKZInfo<Base<F>>
RecursiveBKZ
( Matrix<F>& B,
  Matrix<F>& R,
  Int cutoff,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveBKZ"))
    // TODO: Make this runtime-tunable
    Int numShuffles = 1;
    return bkz::RecursiveHelper( B, R, numShuffles, cutoff, ctrl );
}

template<typename F>
BKZInfo<Base<F>>
BKZ
( Matrix<F>& B,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZ"))
    Matrix<F> R;
    return BKZ( B, R, ctrl );
}

#define PROTO(F) \
  template BKZInfo<Base<F>> BKZ \
  ( Matrix<F>& B, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> BKZ \
  ( Matrix<F>& B, \
    Matrix<F>& R, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> BKZ \
  ( Matrix<F>& B, \
    Matrix<F>& U, \
    Matrix<F>& UInv, \
    Matrix<F>& R, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> BKZWithQ \
  ( Matrix<F>& B, \
    Matrix<F>& QR, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> BKZWithQ \
  ( Matrix<F>& B, \
    Matrix<F>& U, \
    Matrix<F>& UInv, \
    Matrix<F>& QR, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> RecursiveBKZ \
  ( Matrix<F>& B, \
    Int cutoff, \
    const BKZCtrl<Base<F>>& ctrl ); \
  template BKZInfo<Base<F>> RecursiveBKZ \
  ( Matrix<F>& B, \
    Matrix<F>& R, \
    Int cutoff, \
    const BKZCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO // until we have complex enumeration support
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
