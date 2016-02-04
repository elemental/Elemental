/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_BKZ_HPP
#define EL_LATTICE_BKZ_HPP

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
    const Int n = B.Width();

    const Real BOneNorm = OneNorm(B);
    const Real fudge = 2; // TODO: Make tunable
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge));

    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to float");
        try
        {
            typedef float RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, ULower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( U, ULower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to double");
        try
        {
            typedef double RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, ULower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( U, ULower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to DoubleDouble");
        try
        {
            typedef DoubleDouble RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, ULower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( U, ULower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to QuadDouble");
        try
        {
            typedef QuadDouble RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, ULower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( U, ULower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to Quad");
        try
        {
            typedef Quad RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, ULower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( U, ULower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_MPC
    if( !IsFixedPrecision<Real>::value )
    {
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpc::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            if( ctrl.progress )
                Output("Dropping to precision=",neededPrec);
            mpc::SetPrecision( neededPrec );
            try
            {
                Matrix<F> BLower, ULower, QRLower, tLower;
                Matrix<Real> dLower;
                BKZCtrl<Real> ctrlLower( ctrl );
                Copy( B, BLower );
                Copy( U, ULower );
                Copy( QR, QRLower );
                Copy( t, tLower );
                Copy( d, dLower );
                auto infoLower =
                  BKZWithQ
                  ( BLower, ULower, QRLower, tLower, dLower, ctrlLower );
                mpc::SetPrecision( inputPrec );
                BKZInfo<Real> info( infoLower );
                Copy( BLower, B );
                Copy( ULower, U );
                Copy( QRLower, QR );
                Copy( tLower, t );
                Copy( dLower, d );
                return info;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
            mpc::SetPrecision( inputPrec );
        }
    }
#endif

    if( ctrl.recursive &&
        Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n &&
        !ctrl.jumpstart )
    {
        //return RecursiveBKZWithQ( B, U, QR, t, d, ctrl );
        Output("Warning: Computation of U not yet supported for recursive BKZ");
    }

    // TODO: Add optional logging

    if( ctrl.blocksize < 2 )
        LogicError("BKZ requires a blocksize of at least 2");
    if( ctrl.blocksize == 2 )
    {
        auto lllCtrl( ctrl.lllCtrl );
        if( ctrl.jumpstart )
        {
            lllCtrl.jumpstart = true;
            lllCtrl.startCol = 0;
        }
        auto lllInfo = LLLWithQ( B, U, QR, t, d, lllCtrl );
        BKZInfo<Real> info;
        info.delta = lllInfo.delta;
        info.eta = lllInfo.eta;
        info.rank = lllInfo.rank;
        info.nullity = lllInfo.nullity;
        info.numSwaps = lllInfo.numSwaps;
        info.numEnums = 0;
        info.numEnumFailures = 0;
        info.logVol = lllInfo.logVol;
        return info;
    }

    Int numSwaps=0;
    LLLInfo<Real> lllInfo;
    // TODO: Decide on the best way to handle this in the presense of recursion
    //if( ctrl.skipInitialLLL )
    if( false )
    {
        Identity( U, n, n );
        QR = B;
        El::QR( QR, t, d );
    }
    else
    {
        auto lllCtrl( ctrl.lllCtrl );
        if( ctrl.jumpstart )
        {
            lllCtrl.jumpstart = true;
            lllCtrl.startCol = 0;
        }
        lllInfo = LLLWithQ( B, U, QR, t, d, lllCtrl );
        if( ctrl.progress )
            Output("Initial LLL applied ",lllInfo.numSwaps," swaps");
        numSwaps = lllInfo.numSwaps;
    }
    // The zero columns should be at the end of B
    const Int rank = lllInfo.rank;

    ofstream failedEnumFile, streakSizesFile, nontrivialCoordsFile;
    if( ctrl.logFailedEnums )
    {
        Output("Opening failedEnumFile");
        failedEnumFile.open( ctrl.failedEnumFile.c_str() );
    }
    if( ctrl.logStreakSizes )
        streakSizesFile.open( ctrl.streakSizesFile.c_str() );
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.open( ctrl.nontrivialCoordsFile.c_str() );

    Int z=0, j=-1; 
    Int numEnums=0, numEnumFailures=0;
    const Int indent = PushIndent(); 
    while( z < rank-1 ) 
    {
        j = Mod(j+1,rank);
        const Int k = Min(j+ctrl.blocksize-1,rank-1);
        const Int h = Min(k+1,rank-1); 
        
        Matrix<F> v;
        auto BEnum = B( ALL, IR(j,k+1) );
        auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
        const Real oldProjNorm = QR.Get(j,j);
        const Real minProjNorm =
          ShortestVectorEnumeration( BEnum, QREnum, v, ctrl.probabalistic );
        ++numEnums;

        const bool keepMin =
          ( Sqrt(ctrl.lllCtrl.delta)*oldProjNorm > minProjNorm );

        if( keepMin )
        {
            if( ctrl.progress )
            {
                Output
                ("Nontrivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
                Print( v, "v" );
                Output("oldProjNorm=",oldProjNorm,", minProjNorm=",minProjNorm);
            }

            ++numEnumFailures;
            if( ctrl.logFailedEnums )
                failedEnumFile << j << endl;
            if( ctrl.logStreakSizes )
                streakSizesFile << z << endl;
            if( ctrl.logNontrivialCoords )
            {
                for( Int e=0; e<k+1-j; ++e )
                    nontrivialCoordsFile << v.Get(e,0) << " ";
                nontrivialCoordsFile << endl;
            }
            z = 0;

            // Find a unimodular matrix W such that v^T W = [1,0,...,0]
            // and then invert it
            Matrix<F> vTrans, W, Rv;
            Transpose( v, vTrans );
            LLL( vTrans, W, Rv );
            if( vTrans.Get(0,0) == F(1) )
            {
                // Do nothing 
            }
            else if( vTrans.Get(0,0) == F(-1) )
            {
                auto w0 = W( ALL, IR(0) );
                w0 *= Real(-1);
            }
            else
            {
                Print( v, "v" );
                Print( vTrans, "vTrans" );
                Print( W, "W" );
                LogicError("Invalid result of LLL on enumeration coefficients");
            }
            Matrix<F> WInv( W );
            Inverse( WInv );
            Round( WInv );
            // Ensure that we have computed the exact inverse
            Matrix<F> WProd;
            Identity( WProd, k+1-j, k+1-j );
            Gemm( NORMAL, NORMAL, F(-1), W, WInv, F(1), WProd );
            const Real WErr = FrobeniusNorm( WProd );
            if( WErr != Real(0) )
            {
                Print( W, "W" );
                Print( WInv, "invW" );
                LogicError("Did not compute exact inverse of W");
            }

            auto USub = U(ALL,IR(j,k+1));
            auto BSub = B(ALL,IR(j,k+1));
            auto BSubCopy( BSub );
            auto USubCopy( USub );
            Gemm( NORMAL, TRANSPOSE, F(1), BSubCopy, WInv, BSub );
            Gemm( NORMAL, TRANSPOSE, F(1), USubCopy, WInv, USub );
        }
        else
        {
            if( ctrl.progress )
                Output
                ("Trivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
        }

        Matrix<F> W;
        Identity( W, h+1, h+1 );
        
        bool changed = false;
        const auto subInd = IR(0,h+1);
        auto BSub = B( ALL, subInd );
        auto QRSub = QR( ALL, subInd );
        auto tSub = t( subInd, ALL );
        auto dSub = d( subInd, ALL );
        if( ctrl.subBKZ )
        {
            BKZCtrl<Real> subCtrl( ctrl );
            subCtrl.time = false;
            subCtrl.progress = false;
            subCtrl.jumpstart = true;
            // Only if we insist on only one level of recursion
            subCtrl.subBKZ = false;
            subCtrl.blocksize = ctrl.subBlocksizeFunc(ctrl.blocksize);
            subCtrl.earlyAbort = ctrl.subEarlyAbort;
            subCtrl.numEnumsBeforeAbort = ctrl.subNumEnumsBeforeAbort;
            subCtrl.recursive = false;
            subCtrl.logFailedEnums = false;
            subCtrl.logStreakSizes = false;
            subCtrl.logNontrivialCoords = false;
            subCtrl.lllCtrl.jumpstart = false;
            subCtrl.lllCtrl.recursive = false;
            if( ctrl.progress )
              Output("Running sub-BKZ with blocksize=",subCtrl.blocksize);
            auto bkzInfo =
              BKZWithQ( BSub, W, QRSub, tSub, dSub, subCtrl );
            if( ctrl.progress )
              Output
              ("  ",bkzInfo.numSwaps," swaps and ",
               bkzInfo.numEnumFailures," failed enums");
            if( bkzInfo.numSwaps != 0 || bkzInfo.numEnumFailures != 0 )
                changed = true;
            numSwaps += bkzInfo.numSwaps;
        }
        else
        {
            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = ( keepMin ? j : h-1 );
            subLLLCtrl.recursive = false;
            lllInfo = LLLWithQ( BSub, W, QRSub, tSub, dSub, subLLLCtrl );
            if( lllInfo.numSwaps != 0 )
                changed = true;
            numSwaps += lllInfo.numSwaps;
        }
        if( !keepMin )
        {
            if( changed )
            {
                if( ctrl.progress )
                    Output("  Subproblem changed");
                // TODO: Output this j into a file
                z = 0;
            }
            else
                ++z;
        }
        auto USub = U( ALL, subInd );
        auto USubCopy( USub );
        Gemm( NORMAL, NORMAL, F(1), USubCopy, W, USub );

        if( ctrl.earlyAbort &&
            numEnums >= ctrl.numEnumsBeforeAbort &&
            j == rank-1 )
            break;
    }
    SetIndent( indent );

    // Perform a final pass to get the full LLL info
    // NOTE: This could be replaced in favor of manually computing the 
    //       returned LLL info but should be cheap relative to the above BKZ
    LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
    subLLLCtrl.jumpstart = true;
    subLLLCtrl.startCol = n-1;
    subLLLCtrl.recursive = false;
    lllInfo = LLLWithQ( B, U, QR, t, d, subLLLCtrl );
    if( lllInfo.numSwaps != 0 )
        LogicError("Final LLL performed ",lllInfo.numSwaps," swaps");
    numSwaps += lllInfo.numSwaps;

    if( ctrl.logFailedEnums )
        failedEnumFile.close();
    if( ctrl.logStreakSizes )
        streakSizesFile.close();
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.close();

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
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("BKZWithQ"))
    typedef Base<F> Real;
    const Int n = B.Width();

    const Real BOneNorm = OneNorm(B);
    const Real fudge = 2; // TODO: Make tunable
    const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge));

    if( MantissaIsLonger<Real,float>::value &&
        MantissaBits<float>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to float");
        try
        {
            typedef float RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,double>::value &&
        MantissaBits<double>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to double");
        try
        {
            typedef double RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#ifdef EL_HAVE_QD
    if( MantissaIsLonger<Real,DoubleDouble>::value &&
        MantissaBits<DoubleDouble>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to DoubleDouble");
        try
        {
            typedef DoubleDouble RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
    if( MantissaIsLonger<Real,QuadDouble>::value &&
        MantissaBits<QuadDouble>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to QuadDouble");
        try
        {
            typedef QuadDouble RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_QUAD
    if( MantissaIsLonger<Real,Quad>::value &&
        MantissaBits<Quad>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to Quad");
        try
        {
            typedef Quad RealLower;
            typedef ConvertBase<F,RealLower> FLower;
            Matrix<FLower> BLower, QRLower, tLower;
            Matrix<RealLower> dLower;
            BKZCtrl<RealLower> ctrlLower( ctrl );
            Copy( B, BLower );
            Copy( QR, QRLower );
            Copy( t, tLower );
            Copy( d, dLower );
            auto infoLower =
              BKZWithQ( BLower, QRLower, tLower, dLower, ctrlLower );
            BKZInfo<Real> info( infoLower );
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            return info;
        }
        catch( std::exception& e )
        { Output("e.what()=",e.what()); }
    }
#endif
#ifdef EL_HAVE_MPC
    if( !IsFixedPrecision<Real>::value )
    {
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpc::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            if( ctrl.progress )
                Output("Dropping to precision=",neededPrec);
            mpc::SetPrecision( neededPrec );
            try
            {
                Matrix<F> BLower, QRLower, tLower;
                Matrix<Real> dLower;
                BKZCtrl<Real> ctrlLower( ctrl );
                Copy( B, BLower );
                Copy( QR, QRLower );
                Copy( t, tLower );
                Copy( d, dLower );
                auto infoLower =
                  BKZWithQ
                  ( BLower, QRLower, tLower, dLower, ctrlLower );
                mpc::SetPrecision( inputPrec );
                BKZInfo<Real> info( infoLower );
                Copy( BLower, B );
                Copy( QRLower, QR );
                Copy( tLower, t );
                Copy( dLower, d );
                return info;
            }
            catch( std::exception& e )
            { Output("e.what()=",e.what()); }
            mpc::SetPrecision( inputPrec );
        }
    }
#endif

    if( ctrl.recursive &&
        Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n &&
        !ctrl.jumpstart )
        return RecursiveBKZWithQ( B, QR, t, d, ctrl );

    // TODO: Add optional logging

    if( ctrl.blocksize < 2 )
        LogicError("BKZ requires a blocksize of at least 2");
    if( ctrl.blocksize == 2 )
    {
        auto lllCtrl( ctrl.lllCtrl );
        if( ctrl.jumpstart )
        {
            lllCtrl.jumpstart = true;
            lllCtrl.startCol = 0;
        }
        auto lllInfo = LLLWithQ( B, QR, t, d, lllCtrl );
        BKZInfo<Real> info;
        info.delta = lllInfo.delta;
        info.eta = lllInfo.eta;
        info.rank = lllInfo.rank;
        info.nullity = lllInfo.nullity;
        info.numSwaps = lllInfo.numSwaps;
        info.numEnums = 0;
        info.numEnumFailures = 0;
        info.logVol = lllInfo.logVol;
        return info;
    }

    Int numSwaps=0;
    LLLInfo<Real> lllInfo;
    // TODO: Decide on the best way to handle this in the presense of recursion
    //if( ctrl.skipInitialLLL )
    if( false )
    {
        QR = B;
        El::QR( QR, t, d );
    }
    else
    {
        auto lllCtrl( ctrl.lllCtrl );
        if( ctrl.jumpstart )
        {
            lllCtrl.jumpstart = true;
            lllCtrl.startCol = 0;
        }
        lllInfo = LLLWithQ( B, QR, t, d, lllCtrl );
        if( ctrl.progress )
            Output("Initial LLL applied ",lllInfo.numSwaps," swaps");
        numSwaps = lllInfo.numSwaps;
    }
    // The zero columns should be at the end of B
    const Int rank = lllInfo.rank;

    ofstream failedEnumFile, streakSizesFile, nontrivialCoordsFile;
    if( ctrl.logFailedEnums )
    {
        Output("Opening failedEnumFile");
        failedEnumFile.open( ctrl.failedEnumFile.c_str() );
    }
    if( ctrl.logStreakSizes )
        streakSizesFile.open( ctrl.streakSizesFile.c_str() );
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.open( ctrl.nontrivialCoordsFile.c_str() );

    Int z=0, j=-1; 
    Int numEnums=0, numEnumFailures=0;
    const Int indent = PushIndent(); 
    while( z < rank-1 ) 
    {
        j = Mod(j+1,rank);
        const Int k = Min(j+ctrl.blocksize-1,rank-1);
        const Int h = Min(k+1,rank-1); 
        
        Matrix<F> v;
        auto BEnum = B( ALL, IR(j,k+1) );
        auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
        const Real oldProjNorm = QR.Get(j,j);
        const Real minProjNorm =
          ShortestVectorEnumeration( BEnum, QREnum, v, ctrl.probabalistic );
        ++numEnums;

        const bool keepMin =
          ( Sqrt(ctrl.lllCtrl.delta)*oldProjNorm > minProjNorm );

        if( keepMin )
        {
            if( ctrl.progress )
            {
                Output
                ("Nontrivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
                Print( v, "v" );
                Output("oldProjNorm=",oldProjNorm,", minProjNorm=",minProjNorm);
            }

            ++numEnumFailures;
            if( ctrl.logFailedEnums )
                failedEnumFile << j << endl;
            if( ctrl.logStreakSizes )
                streakSizesFile << z << endl;
            if( ctrl.logNontrivialCoords )
            {
                for( Int e=0; e<k+1-j; ++e )
                    nontrivialCoordsFile << v.Get(e,0) << " ";
                nontrivialCoordsFile << endl;
            }
            z = 0;

            // Find a unimodular matrix W such that v^T W = [1,0,...,0]
            // and then invert it
            Matrix<F> vTrans, W, Rv;
            Transpose( v, vTrans );
            LLL( vTrans, W, Rv );
            if( vTrans.Get(0,0) == F(1) )
            {
                // Do nothing 
            }
            else if( vTrans.Get(0,0) == F(-1) )
            {
                auto w0 = W( ALL, IR(0) );
                w0 *= Real(-1);
            }
            else
            {
                Print( v, "v" );
                Print( vTrans, "vTrans" );
                Print( W, "W" );
                LogicError("Invalid result of LLL on enumeration coefficients");
            }
            Matrix<F> WInv( W );
            Inverse( WInv );
            Round( WInv );
            // Ensure that we have computed the exact inverse
            Matrix<F> WProd;
            Identity( WProd, k+1-j, k+1-j );
            Gemm( NORMAL, NORMAL, F(-1), W, WInv, F(1), WProd );
            const Real WErr = FrobeniusNorm( WProd );
            if( WErr != Real(0) )
            {
                Print( W, "W" );
                Print( WInv, "invW" );
                LogicError("Did not compute exact inverse of W");
            }

            auto BSub = B(ALL,IR(j,k+1));
            auto BSubCopy( BSub );
            Gemm( NORMAL, TRANSPOSE, F(1), BSubCopy, WInv, BSub );
        }
        else
        {
            if( ctrl.progress )
                Output
                ("Trivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
        }
        
        bool changed = false;
        const auto subInd = IR(0,h+1);
        auto BSub = B( ALL, subInd );
        auto QRSub = QR( ALL, subInd );
        auto tSub = t( subInd, ALL );
        auto dSub = d( subInd, ALL );
        if( ctrl.subBKZ )
        {
            BKZCtrl<Real> subCtrl( ctrl );
            subCtrl.time = false;
            subCtrl.progress = false;
            subCtrl.jumpstart = true;
            // Only if we insist on only one level of recursion
            subCtrl.subBKZ = false;
            subCtrl.blocksize = ctrl.subBlocksizeFunc(ctrl.blocksize);
            subCtrl.earlyAbort = ctrl.subEarlyAbort;
            subCtrl.numEnumsBeforeAbort = ctrl.subNumEnumsBeforeAbort;
            subCtrl.recursive = false;
            subCtrl.logFailedEnums = false;
            subCtrl.logStreakSizes = false;
            subCtrl.logNontrivialCoords = false;
            subCtrl.lllCtrl.jumpstart = false;
            subCtrl.lllCtrl.recursive = false;
            if( ctrl.progress )
              Output("Running sub-BKZ with blocksize=",subCtrl.blocksize);
            auto bkzInfo = BKZWithQ( BSub, QRSub, tSub, dSub, subCtrl );
            if( ctrl.progress )
              Output
              ("  ",bkzInfo.numSwaps," swaps and ",
               bkzInfo.numEnumFailures," failed enums");
            if( bkzInfo.numSwaps != 0 || bkzInfo.numEnumFailures != 0 )
                changed = true;
            numSwaps += bkzInfo.numSwaps;
        }
        else
        {
            LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
            subLLLCtrl.jumpstart = true;
            subLLLCtrl.startCol = ( keepMin ? j : h-1 );
            subLLLCtrl.recursive = false;
            lllInfo = LLLWithQ( BSub, QRSub, tSub, dSub, subLLLCtrl );
            if( lllInfo.numSwaps != 0 )
                changed = true;
            numSwaps += lllInfo.numSwaps;
        }
        if( !keepMin )
        {
            if( changed )
            {
                if( ctrl.progress )
                    Output("  Subproblem changed");
                // TODO: Output this j into a file
                z = 0;
            }
            else
                ++z;
        }

        if( ctrl.earlyAbort &&
            numEnums >= ctrl.numEnumsBeforeAbort &&
            j == rank-1 )
            break;
    }
    SetIndent( indent );

    // Perform a final pass to get the full LLL info
    // NOTE: This could be replaced in favor of manually computing the 
    //       returned LLL info but should be cheap relative to the above BKZ
    LLLCtrl<Real> subLLLCtrl( ctrl.lllCtrl );
    subLLLCtrl.jumpstart = true;
    subLLLCtrl.startCol = n-1;
    subLLLCtrl.recursive = false;
    lllInfo = LLLWithQ( B, QR, t, d, subLLLCtrl );
    if( lllInfo.numSwaps != 0 )
        LogicError("Final LLL performed ",lllInfo.numSwaps," swaps");
    numSwaps += lllInfo.numSwaps;

    if( ctrl.logFailedEnums )
        failedEnumFile.close();
    if( ctrl.logStreakSizes )
        streakSizesFile.close();
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.close();

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
        const Real fudge = 2; // TODO: Make tunable
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
        if( !succeeded && !IsFixedPrecision<Real>::value )
        {
            // Only move down to a lower-precision MPFR type if the jump is
            // substantial. The current value has been naively chosen.
            const mpfr_prec_t minPrecDiff = 32;
            mpfr_prec_t inputPrec = mpc::Precision();
            if( neededPrec <= inputPrec-minPrecDiff )
            {
                mpc::SetPrecision( neededPrec );
                try {
                    info =
                      LowerPrecisionMerge<F,BigFloat>
                      ( CL, CR, B, QR, t, d, ctrl );
                    info.numSwaps += numPrevSwaps;
                    succeeded = true;
                }
                catch( std::exception& e )
                { Output("e.what()=",e.what()); }
                mpc::SetPrecision( inputPrec );
            }
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
        const Real fudge = 2; // TODO: Make tunable
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
#ifdef EL_HAVE_MPC
        if( !succeeded && !IsFixedPrecision<Real>::value )
        {
            // Only move down to a lower-precision MPFR type if the jump is
            // substantial. The current value has been naively chosen.
            const mpfr_prec_t minPrecDiff = 32;
            mpfr_prec_t inputPrec = mpc::Precision();
            if( neededPrec <= inputPrec-minPrecDiff )
            {
                mpc::SetPrecision( neededPrec );
                try {
                    info =
                      LowerPrecisionMerge<F,BigFloat>
                      ( CL, CR, B, QR, t, d, ctrl );
                    info.numSwaps += numPrevSwaps;
                    succeeded = true;
                }
                catch( std::exception& e )
                { Output("e.what()=",e.what()); }
                mpc::SetPrecision( inputPrec );
            }
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

} // namespace El

#endif // ifndef EL_LATTICE_BKZ_HPP
