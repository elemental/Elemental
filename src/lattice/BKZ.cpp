/*
   Copyright (c) 2009-2016, Jack Poulson
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
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZWithQ"))
    typedef Base<F> Real;
    const Int m = B.Height();
    const Int n = B.Width();
    if( ctrl.recursive && Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n )
    {
        //return RecursiveBKZWithQ( B, U, QR, t, d, ctrl );
        Output("Warning: Computation of U not yet supported for recursive BKZ");
    }

    if( ctrl.blocksize < 2 )
        LogicError("BKZ requires a blocksize of at least 2");
    if( m < n )
        LogicError("If height(B) < width(B), call LLL first");

    const bool progress = true;

    Int numSwaps=0;
    auto lllInfo = LLLWithQ( B, U, QR, t, d, ctrl.lllCtrl );
    numSwaps = lllInfo.numSwaps;

    Int z=0, j=-1; 
    Matrix<F> BTmp, UTmp, QRTmp, tTmp;
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
            auto QRSub = QR( ALL, subInd );
            auto tSub = t( subInd, ALL );
            auto dSub = d( subInd, ALL );
            lllInfo = LLLWithQ( BSub, USub, QRSub, tSub, dSub, subLLLCtrl );
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
            // The initializations of the expanded U is simpler to derive.
 
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
            lllInfo = LLLWithQ( BTmp, UTmp, QRTmp, tTmp, dTmp, subLLLCtrl );
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
            }
            // Returning the QR factorization doesn't work without explicitly
            // forming Q due to the first column of BTmp being removed
            subLLLCtrl.startCol = 0;
            lllInfo = LLLWithQ( B, U, QR, t, d, ctrl.lllCtrl );
        }

        if( ctrl.earlyAbort && numEnums >= ctrl.numEnumsBeforeAbort )
            break;
    }

    // Perform a final pass to get the full LLL info
    // NOTE: This could be replaced in favor of manually computing the 
    //       returned LLL info but should be cheap relative to the above BKZ
    lllInfo = LLLWithQ( B, U, QR, t, d, ctrl.lllCtrl );
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
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZ"))
    typedef Base<F> Real;
    Matrix<F> t;
    Matrix<Real> d;
    auto info = BKZWithQ( B, U, R, t, d, ctrl );
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
    if( ctrl.recursive && Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n )
        return RecursiveBKZWithQ( B, QR, t, d, ctrl );

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
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
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
        auto cL = CL( ALL, IR(jSub) );
        auto cR = CR( ALL, IR(jSub) );
        auto bL = BLower( ALL, IR(2*jSub) );
        auto bR = BLower( ALL, IR(2*jSub+1) );
        Copy( cL, bL );
        Copy( cR, bR );
    }
    if( firstHalf > n/2 )
    {
        auto cL = CL( ALL, IR(firstHalf-1) );
        auto bL = BLower( ALL, IR(n-1) );
        Copy( cL, bL );
    }

    BKZCtrl<RealLower> ctrlLower( ctrl );
    ctrlLower.recursive = false;
    ctrlLower.lllCtrl.recursive = false;
    RealLower eps = limits::Epsilon<RealLower>();
    RealLower minEta = RealLower(1)/RealLower(2)+Pow(eps,RealLower(0.9));
    ctrlLower.lllCtrl.eta = Max(minEta,ctrlLower.lllCtrl.eta);
    Timer timer;
    Matrix<FLower> QRLower, tLower;
    Matrix<RealLower> dLower;
    if( ctrl.lllCtrl.time )
        timer.Start();
    auto infoLower = BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
    if( ctrl.lllCtrl.time )
        Output("  " + typeString + " BKZ took ",timer.Stop()," seconds");
    Copy( BLower, B );
    Copy( QRLower, QR );
    Copy( tLower, t ); 
    Copy( dLower, d );
    return infoLower;
}

template<typename Real>
BKZInfo<Real>
RecursiveHelper
( Matrix<Real>& B,
  Matrix<Real>& U,
  Matrix<Real>& QR,
  Matrix<Real>& t,
  Matrix<Real>& d,
  Int numShuffles,
  bool maintainU,
  const BKZCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bkz::RecursiveHelper"))
    if( maintainU )
        LogicError("Recursive BKZ does not yet support computing U");

    typedef Real F;
    const Int n = B.Width();
    if( n <= Max(ctrl.lllCtrl.cutoff,ctrl.blocksize) )
    {
        auto ctrlMod( ctrl );
        ctrlMod.recursive = false;
        ctrlMod.lllCtrl.recursive = false;
        if( maintainU )
            return BKZWithQ( B, U, QR, t, d, ctrlMod );
        else
            return BKZWithQ( B, QR, t, d, ctrlMod );
    }
    Timer timer;

    // Reduce the entire matrix with LLL before attempting BKZ.
    // Deep reductions should probably not be used due to the expense.
    BKZInfo<Real> info;
    info.numSwaps = 0;
    auto lllCtrlMod( ctrl.lllCtrl );
    lllCtrlMod.recursive = true;
    auto lllInfo = LLL( B, QR, lllCtrlMod );
    info.numSwaps += lllInfo.numSwaps;

    const Real BOneNorm = OneNorm( B );
    Output("|| B ||_1 = ",BOneNorm);

    for( Int shuffle=0; shuffle<=numShuffles; ++shuffle )
    {
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
            Output("Shuffle=",shuffle);
        auto C( B ); 

        const Int firstHalf = n-(n/2);
        Range<Int> indL(0,firstHalf), indR(firstHalf,n);
        auto CL = C( ALL, indL );
        auto CR = C( ALL, indR );

        double leftTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        BKZInfo<Real> leftInfo;
        {
            Matrix<Real> QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveBKZWithQ( CL, QRL, tL, dL, ctrl ); 
        }
        if( ctrl.lllCtrl.time )
            leftTime = timer.Stop(); 

        double rightTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        BKZInfo<Real> rightInfo;
        {
            Matrix<Real> QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveBKZWithQ( CR, QRR, tR, dR, ctrl );
        }
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
        const unsigned neededPrec = unsigned(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( MantissaIsLonger<Real,float>::value &&
            MantissaBits<float>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,float>( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            MantissaIsLonger<Real,double>::value &&
            MantissaBits<double>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,double>( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#ifdef EL_HAVE_QD
        if( !succeeded &&
            MantissaIsLonger<Real,DoubleDouble>::value &&
            MantissaBits<DoubleDouble>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,DoubleDouble>
                  ( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
        if( !succeeded &&
            MantissaIsLonger<Real,QuadDouble>::value &&
            MantissaBits<QuadDouble>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,QuadDouble>
                  ( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
        }
#elif defined(EL_HAVE_QUAD)
        if( !succeeded &&
            MantissaIsLonger<Real,Quad>::value &&
            MantissaBits<Quad>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,Quad>( CL, CR, B, QR, t, d, ctrl );
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
                info =
                  LowerPrecisionMerge<F,BigFloat>( CL, CR, B, QR, t, d, ctrl );
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
            
            auto ctrlMod( ctrl );
            ctrlMod.recursive = false;
            ctrlMod.lllCtrl.recursive = false;
            info = BKZWithQ( B, QR, t, d, ctrlMod );
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
  Matrix<Complex<Real>>& U,
  Matrix<Complex<Real>>& QR,
  Matrix<Complex<Real>>& t,
  Matrix<Real>& d,
  Int numShuffles,
  bool maintainU,
  const BKZCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("bkz::RecursiveHelper"))
    if( maintainU )
        LogicError("Recursive BKZ does not yet support computing U");

    typedef Complex<Real> F;
    const Int n = B.Width();
    if( n <= Max(ctrl.lllCtrl.cutoff,ctrl.blocksize) )
    {
        auto ctrlMod( ctrl );
        ctrlMod.recursive = false;
        ctrlMod.lllCtrl.recursive = false;
        if( maintainU )
            return BKZWithQ( B, U, QR, t, d, ctrlMod );
        else
            return BKZWithQ( B, QR, t, d, ctrlMod );
    }
    Timer timer;

    // Reduce the entire matrix with LLL before attempting BKZ.
    // Deep reductions should probably not be used due to the expense.
    BKZInfo<Real> info;
    info.numSwaps = 0;
    auto lllCtrlMod( ctrl.lllCtrl );
    lllCtrlMod.recursive = true;
    auto lllInfo = LLL( B, QR, lllCtrlMod );
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

        double leftTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        LLLInfo<Real> leftInfo;
        {
            Matrix<Complex<Real>> QRL, tL;
            Matrix<Real> dL;
            leftInfo = RecursiveBKZWithQ( CL, QRL, tL, dL, ctrl ); 
        }
        if( ctrl.lllCtrl.time )
            leftTime = timer.Stop(); 

        double rightTime;
        if( ctrl.lllCtrl.time )
            timer.Start();
        LLLInfo<Real> rightInfo;
        {
            Matrix<Complex<Real>> QRR, tR;
            Matrix<Real> dR;
            rightInfo = RecursiveBKZWithQ( CR, QRR, tR, dR, ctrl );
        }
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
        const unsigned neededPrec = unsigned(Ceil(Log2(COneNorm)*fudge));
        if( ctrl.lllCtrl.progress || ctrl.lllCtrl.time )
        {
            Output("  || C ||_1 = ",COneNorm);
            Output("  Needed precision: ",neededPrec);
        }

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        if( MantissaIsLonger<Real,float>::value &&
            MantissaBits<float>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,float>( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
        if( !succeeded && 
            MantissaIsLonger<Real,double>::value &&
            MantissaBits<double>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,double>( CL, CR, B, QR, t, d, ctrl );
                info.numSwaps += numPrevSwaps;
                succeeded = true;
            } catch( std::exception& e ) { Output("e.what()=",e.what()); }
        }
        // There is not yet support for Complex<{Quad,Double}Double>
#ifdef EL_HAVE_QUAD
        if( !succeeded &&
            MantissaIsLonger<Real,Quad>::value &&
            MantissaBits<Quad>::value >= neededPrec )
        {
            try
            {
                info =
                  LowerPrecisionMerge<F,Quad>( CL, CR, B, QR, t, d, ctrl );
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
            
            auto ctrlMod( ctrl );
            ctrlMod.recursive = false;
            ctrlMod.lllCtrl.recursive = false;
            info = BKZWithQ( B, QR, t, d, ctrlMod );
            info.numSwaps += numPrevSwaps;
        }
    }
    return info;
}

} // namespace bkz

template<typename F>
BKZInfo<Base<F>>
RecursiveBKZWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveBKZWithQ"))
    // TODO: Make this runtime-tunable
    Int numShuffles = 1;
    Matrix<F> U;
    bool maintainU=false;
    return bkz::RecursiveHelper( B, U, QR, t, d, numShuffles, maintainU, ctrl );
}

template<typename F>
BKZInfo<Base<F>>
RecursiveBKZWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RecursiveBKZWithQ"))
    // TODO: Make this runtime-tunable
    Int numShuffles = 1;
    bool maintainU=true;
    return bkz::RecursiveHelper( B, U, QR, t, d, numShuffles, maintainU, ctrl );
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
    Matrix<F>& QR, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    const BKZCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO // until we have complex enumeration support
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
