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

// Block Korkin-Zolotarev (BKZ) reduction
// ======================================
// TODO: Tailor BKZInfo; it is currently a copy of LLLInfo
template<typename Real>
struct BKZInfo
{
    Real delta;
    Real eta; 
    Int rank;
    Int nullity; 
    Int numSwaps;
    Int numEnums;
    Int numEnumFailures;
    Real logVol;

    template<typename OtherReal>
    BKZInfo<Real>& operator=( const BKZInfo<OtherReal>& info )
    {
        delta = Real(info.delta);
        eta = Real(info.eta);
        rank = info.rank;
        nullity = info.nullity;
        numSwaps = info.numSwaps;
        numEnums = info.numEnums;
        numEnumFailures = info.numEnumFailures;
        logVol = Real(info.logVol);
        return *this;
    }

    BKZInfo() { }
    BKZInfo( const BKZInfo<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    BKZInfo( const BKZInfo<OtherReal>& ctrl ) { *this = ctrl; }
};

template<typename Real>
struct BKZCtrl
{
    Int blocksize=20;
    bool time=false;
    bool progress=false;

    bool earlyAbort=false;
    Int numEnumsBeforeAbort=1000; // only used if earlyAbort=true

    bool variableBlocksize=false;
    function<Int(Int)> blocksizeFunc;

    bool variableEnumType=false;
    function<EnumType(Int)> enumTypeFunc;

    // Y-sparse enumeration supports simultaneous searches for improving a
    // contiguous window of vectors
    Int multiEnumWindow=15;

    bool skipInitialLLL=false;
    bool jumpstart=false;
    Int startCol=0;

    EnumCtrl<Real> enumCtrl;

    // Rather than running LLL after a productive enumeration, one could run
    // BKZ with a smaller blocksize (perhaps with early abort)
    bool subBKZ=true;
    function<Int(Int)> subBlocksizeFunc =
      function<Int(Int)>( []( Int bsize )
      { return Max(Min(bsize/2,Int(20)),Int(2)); } );
    bool subEarlyAbort = true;
    Int subNumEnumsBeforeAbort = 100;

    // This seems to be *more* expensive but lead to higher quality (perhaps).
    // Note that this is different from GNR recursion and uses a tree method
    // with shuffling at each merge.
    bool recursive=false;

    bool logFailedEnums=false;
    std::string failedEnumFile="BKZFailedEnums.txt";

    bool logStreakSizes=false;
    std::string streakSizesFile="BKZStreakSizes.txt";

    bool logNontrivialCoords=false;
    std::string nontrivialCoordsFile="BKZNontrivialCoords.txt";

    bool logNorms=false;
    std::string normsFile="BKZNorms.txt";

    bool logProjNorms=false;
    std::string projNormsFile="BKZProjNorms.txt";

    bool checkpoint=false;
    FileFormat checkpointFormat=ASCII;
    std::string checkpointFileBase="BKZCheckpoint";
    std::string tourFileBase="BKZTour";

    LLLCtrl<Real> lllCtrl;

    // We frequently need to convert datatypes, so make this easy
    template<typename OtherReal>
    BKZCtrl<Real>& operator=( const BKZCtrl<OtherReal>& ctrl )
    {
        blocksize = ctrl.blocksize;
        time = ctrl.time;
        progress = ctrl.progress;

        earlyAbort = ctrl.earlyAbort;
        numEnumsBeforeAbort = ctrl.numEnumsBeforeAbort;

        variableBlocksize = ctrl.variableBlocksize;
        blocksizeFunc = ctrl.blocksizeFunc;

        variableEnumType = ctrl.variableEnumType;
        enumTypeFunc = ctrl.enumTypeFunc;

        multiEnumWindow = ctrl.multiEnumWindow;

        skipInitialLLL = ctrl.skipInitialLLL;
        jumpstart = ctrl.jumpstart;
        startCol = ctrl.startCol;

        enumCtrl = ctrl.enumCtrl;

        subBKZ = ctrl.subBKZ;
        subBlocksizeFunc = ctrl.subBlocksizeFunc;
        subEarlyAbort = ctrl.subEarlyAbort;
        subNumEnumsBeforeAbort = ctrl.subNumEnumsBeforeAbort;

        recursive = ctrl.recursive;

        logFailedEnums = ctrl.logFailedEnums;
        logStreakSizes = ctrl.logStreakSizes;
        logNontrivialCoords = ctrl.logNontrivialCoords;
        logNorms = ctrl.logNorms;
        logProjNorms = ctrl.logProjNorms;
        checkpoint = ctrl.checkpoint;
        failedEnumFile = ctrl.failedEnumFile;
        streakSizesFile = ctrl.streakSizesFile;
        nontrivialCoordsFile = ctrl.nontrivialCoordsFile;
        normsFile = ctrl.normsFile;
        projNormsFile = ctrl.projNormsFile;
        checkpointFileBase = ctrl.checkpointFileBase;
        tourFileBase = ctrl.tourFileBase;
        checkpointFormat = ctrl.checkpointFormat;

        lllCtrl = ctrl.lllCtrl;
        return *this;
    }

    BKZCtrl() { }
    BKZCtrl( const BKZCtrl<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    BKZCtrl( const BKZCtrl<OtherReal>& ctrl ) { *this = ctrl; }
};

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& R,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

template<typename F>
BKZInfo<Base<F>> BKZWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl=BKZCtrl<Base<F>>() );

namespace bkz {

static Timer enumTimer, bkzTimer;

template<typename F>
bool TrivialCoordinates( const Matrix<F>& v )
{
    EL_DEBUG_CSE
    const Int n = v.Height();    
    if( n == 0 )
        LogicError("Invalid coordinate length");
    if( v(0) != F(1) )
        return false;
    for( Int i=1; i<n; ++i )
        if( v(i) != F(0) )
            return false;
    return true;
}

template<typename RealLower,typename F>
bool TryLowerPrecision
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
  unsigned neededPrec,
  BKZInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( MantissaIsLonger<Base<F>,RealLower>::value &&
        MantissaBits<RealLower>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to ",TypeName<RealLower>());
        try
        {
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
            info = infoLower;
            Copy( BLower, B );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            succeeded = true;
        }
        catch( std::exception& e )
        {
            Output("e.what()=",e.what());
        }
    }
    return succeeded;
}

template<typename RealLower,typename F>
bool TryLowerPrecision
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
  unsigned neededPrec,
  BKZInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( MantissaIsLonger<Base<F>,RealLower>::value &&
        MantissaBits<RealLower>::value >= neededPrec )
    {
        if( ctrl.progress )
            Output("Dropping to ",TypeName<RealLower>());
        try
        {
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
            info = infoLower;
            Copy( BLower, B );
            Copy( ULower, U );
            Copy( QRLower, QR );
            Copy( tLower, t );
            Copy( dLower, d );
            succeeded = true;
        }
        catch( std::exception& e )
        {
            Output("e.what()=",e.what());
        }
    }
    return succeeded;
}

#ifdef EL_HAVE_MPC
template<typename F>
bool TryLowerPrecisionBigFloat
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
  unsigned neededPrec,
  BKZInfo<Base<F>>& info )
{
    typedef Base<F> Real;
    bool succeeded = false;
    if( !IsFixedPrecision<Real>::value )
    {
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpfr::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            if( ctrl.progress )
                Output("Dropping to precision=",neededPrec);
            mpfr::SetPrecision( neededPrec );
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
                mpfr::SetPrecision( inputPrec );
                info = infoLower;
                Copy( BLower, B );
                Copy( QRLower, QR );
                Copy( tLower, t );
                Copy( dLower, d );
                succeeded = true;
            }
            catch( std::exception& e )
            {
                Output("e.what()=",e.what());
            }
            mpfr::SetPrecision( inputPrec );
        }
    }
    return succeeded;
}

template<typename F>
bool TryLowerPrecisionBigFloat
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
  unsigned neededPrec,
  BKZInfo<Base<F>>& info )
{
    typedef Base<F> Real;
    bool succeeded = false;
    if( !IsFixedPrecision<Real>::value )
    {
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpfr::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            if( ctrl.progress )
                Output("Dropping to precision=",neededPrec);
            mpfr::SetPrecision( neededPrec );
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
                mpfr::SetPrecision( inputPrec );
                info = infoLower;
                Copy( BLower, B );
                Copy( ULower, U );
                Copy( QRLower, QR );
                Copy( tLower, t );
                Copy( dLower, d );
                succeeded = true;
            }
            catch( std::exception& e )
            {
                Output("e.what()=",e.what());
            }
            mpfr::SetPrecision( inputPrec );
        }
    }
    return succeeded;
}
#endif

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
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = B.Width();

    const bool isInteger = IsInteger( B );
    if( isInteger )
    {
        const Real BOneNorm = OneNorm(B);
        const Real fudge = 2; // TODO: Make tunable
        const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge));

        BKZInfo<Real> info;
        bool succeeded = bkz::TryLowerPrecision<float>
          ( B, U, QR, t, d, ctrl, neededPrec, info );
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<double>
              ( B, U, QR, t, d, ctrl, neededPrec, info );
#ifdef EL_HAVE_QD
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<DoubleDouble>
              ( B, U, QR, t, d, ctrl, neededPrec, info );
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<QuadDouble>
              ( B, U, QR, t, d, ctrl, neededPrec, info );
#elif defined(EL_HAVE_QUAD)
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<Quad>
              ( B, U, QR, t, d, ctrl, neededPrec, info );
#endif
#ifdef EL_HAVE_MPC
        if( !succeeded )
            succeeded = bkz::TryLowerPrecisionBigFloat
              ( B, U, QR, t, d, ctrl, neededPrec, info );
#endif
        if( succeeded )
            return info;
    }
    // TODO: Allow for dropping with non-integer vectors?

    if( ctrl.recursive &&
        Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n &&
        !ctrl.jumpstart )
    {
        //return RecursiveBKZWithQ( B, U, QR, t, d, ctrl );
        Output("Warning: Computation of U not yet supported for recursive BKZ");
    }

    if( ctrl.time )
    {
        bkz::enumTimer.Reset();
        bkz::bkzTimer.Reset();
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
        auto lllInfo = LLLWithQ( B, U, QR, t, d, lllCtrl );
        if( ctrl.time )
            Output("Initial LLL time: ",bkz::bkzTimer.Stop()," seconds");
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
        lllInfo = LLLWithQ( B, U, QR, t, d, lllCtrl );
        if( ctrl.progress )
            Output("Initial LLL applied ",lllInfo.numSwaps," swaps");
        if( ctrl.time )
            Output("Initial LLL time: ",bkz::bkzTimer.Stop()," seconds");
        numSwaps = lllInfo.numSwaps;
    }
    // The zero columns should be at the end of B
    const Int rank = lllInfo.rank;

    ofstream failedEnumFile, streakSizesFile,
             normsFile, projNormsFile,
             nontrivialCoordsFile;
    if( ctrl.logFailedEnums )
        failedEnumFile.open( ctrl.failedEnumFile.c_str() );
    if( ctrl.logStreakSizes )
        streakSizesFile.open( ctrl.streakSizesFile.c_str() );
    if( ctrl.logNorms )
        normsFile.open( ctrl.normsFile.c_str() );
    if( ctrl.logProjNorms )
        projNormsFile.open( ctrl.projNormsFile.c_str() );
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.open( ctrl.nontrivialCoordsFile.c_str() );

    auto enumCtrl = ctrl.enumCtrl;
    enumCtrl.disablePrecDrop = true;

    Int z=0;
    Int j = ( ctrl.jumpstart ? ctrl.startCol : 0 ) - 1;
    Int numEnums=0, numEnumFailures=0;
    const Int indent = PushIndent(); 
    while( z < rank-1 ) 
    {
        j = Mod(j+1,rank);
        Int bsize = ctrl.blocksize;
        if( ctrl.variableBlocksize )
            bsize = ctrl.blocksizeFunc(j);
        const Int k = Min(j+bsize-1,rank-1);
        const Int h = Min(k+1,rank-1); 
        if( ctrl.checkpoint )
        {
            Write( B, ctrl.checkpointFileBase, ctrl.checkpointFormat, "B" );
            if( j == 0 )
                Write( B, ctrl.tourFileBase, ctrl.checkpointFormat, "B" );
        }
        if( j == 0 )
        {
            if( ctrl.logNorms )
            {
                for( Int j=0; j<n; ++j )
                    normsFile << FrobeniusNorm(B(ALL,IR(j))) << " ";    
                normsFile << endl;
            }
            if( ctrl.logProjNorms )
            {
                for( Int j=0; j<n; ++j )
                    projNormsFile << RealPart(QR(j,j)) << " ";    
                projNormsFile << endl;
            }
        }
        
        Matrix<F> v;
        auto BEnum = B( ALL, IR(j,k+1) );
        auto UEnum = U( ALL, IR(j,k+1) );
        auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
        if( ctrl.time )
            bkz::enumTimer.Start();
        if( ctrl.variableEnumType )
            enumCtrl.enumType = ctrl.enumTypeFunc(j);
        const Range<Int> windowInd = IR(j,Min(j+ctrl.multiEnumWindow,k+1));
        auto normUpperBounds = GetRealPartOfDiagonal(QR(windowInd,windowInd));
        Scale( Min(Sqrt(ctrl.lllCtrl.delta),Real(1)), normUpperBounds );
        const auto minPair = 
          MultiShortestVectorEnrichment
          ( BEnum, UEnum, QREnum, normUpperBounds, v, enumCtrl );
        if( ctrl.time )
            Output("Enum/enrich time: ",bkz::enumTimer.Stop()," seconds");
        ++numEnums;

        const Real minProjNorm = minPair.first;
        const Int insertionInd = minPair.second;

        const Real oldProjNorm = RealPart(QREnum(insertionInd,insertionInd));
        const bool keptMin = ( minProjNorm < oldProjNorm );
        if( keptMin )
        {
            if( ctrl.progress )
            {
                Output
                ("Nontrivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
                Print( v, "v" );
                Output("insertion index: ",insertionInd);
                Output("oldProjNorm=",oldProjNorm,", minProjNorm=",minProjNorm,", oldProjNorm-minProjNorm=",oldProjNorm-minProjNorm);
            }

            ++numEnumFailures;
            if( ctrl.logFailedEnums )
                failedEnumFile << j << endl;
            if( ctrl.logStreakSizes )
                streakSizesFile << z << endl;
            if( ctrl.logNontrivialCoords )
            {
                for( Int e=0; e<v.Height(); ++e )
                    nontrivialCoordsFile << v(e) << " ";
                nontrivialCoordsFile << endl;
            }
            z = 0;
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
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
            subCtrl.variableBlocksize = false;
            subCtrl.variableEnumType = false;
            subCtrl.recursive = false;
            subCtrl.logFailedEnums = false;
            subCtrl.logStreakSizes = false;
            subCtrl.logNontrivialCoords = false;
            subCtrl.logNorms = false;
            subCtrl.logProjNorms = false;
            subCtrl.checkpoint = false;
            subCtrl.enumCtrl.disablePrecDrop = true;
            subCtrl.enumCtrl.time = false;
            subCtrl.enumCtrl.progress = false;
            subCtrl.lllCtrl.jumpstart = false;
            subCtrl.lllCtrl.recursive = false;
            if( ctrl.progress )
              Output("Running sub-BKZ with blocksize=",subCtrl.blocksize);
            auto bkzInfo = BKZWithQ( BSub, W, QRSub, tSub, dSub, subCtrl );
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
            subLLLCtrl.startCol = ( keptMin ? j : h-1 );
            subLLLCtrl.recursive = false;
            lllInfo = LLLWithQ( BSub, W, QRSub, tSub, dSub, subLLLCtrl );
            if( lllInfo.numSwaps != 0 )
                changed = true;
            numSwaps += lllInfo.numSwaps;
        }
        auto USub = U( ALL, subInd );
        auto USubCopy( USub );
        Gemm( NORMAL, NORMAL, F(1), USubCopy, W, USub );
        if( ctrl.time )
            Output("BKZ time: ",bkz::bkzTimer.Stop()," seconds");
        if( !keptMin )
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
    lllInfo = LLLWithQ( B, U, QR, t, d, subLLLCtrl );
    if( lllInfo.numSwaps != 0 )
        LogicError("Final LLL performed ",lllInfo.numSwaps," swaps");
    numSwaps += lllInfo.numSwaps;

    if( ctrl.logFailedEnums )
        failedEnumFile.close();
    if( ctrl.logStreakSizes )
        streakSizesFile.close();
    if( ctrl.logNorms )
        normsFile.close();
    if( ctrl.logProjNorms )
        projNormsFile.close();
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.close();

    if( ctrl.time )
    {
        Output("Total enumeration time: ",bkz::enumTimer.Total()," seconds");
        Output("Total sub-BKZ time:     ",bkz::bkzTimer.Total()," seconds");
    }

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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = B.Width();

    const bool isInteger = IsInteger(B);
    if( isInteger )
    {
        const Real BOneNorm = OneNorm(B);
        const Real fudge = 2; // TODO: Make tunable
        const unsigned neededPrec = unsigned(Ceil(Log2(BOneNorm)*fudge));

        BKZInfo<Real> info;
        bool succeeded = bkz::TryLowerPrecision<float>
          ( B, QR, t, d, ctrl, neededPrec, info );
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<double>
              ( B, QR, t, d, ctrl, neededPrec, info );
#ifdef EL_HAVE_QD
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<DoubleDouble>
              ( B, QR, t, d, ctrl, neededPrec, info );
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<QuadDouble>
              ( B, QR, t, d, ctrl, neededPrec, info );
#elif defined(EL_HAVE_QUAD)
        if( !succeeded )
            succeeded = bkz::TryLowerPrecision<Quad>
              ( B, QR, t, d, ctrl, neededPrec, info );
#endif
#ifdef EL_HAVE_MPC
        if( !succeeded )
            succeeded = bkz::TryLowerPrecisionBigFloat
              ( B, QR, t, d, ctrl, neededPrec, info );
#endif
        if( succeeded )
            return info;
    }
    // TODO: Allow for dropping with non-integer vectors?

    if( ctrl.recursive &&
        Max(ctrl.blocksize,ctrl.lllCtrl.cutoff) < n &&
        !ctrl.jumpstart )
        return RecursiveBKZWithQ( B, QR, t, d, ctrl );

    if( ctrl.time )
    {
        bkz::enumTimer.Reset();
        bkz::bkzTimer.Reset();
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
        auto lllInfo = LLLWithQ( B, QR, t, d, lllCtrl );
        if( ctrl.time )
            Output("Initial LLL time: ",bkz::bkzTimer.Stop()," seconds");
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
        lllInfo = LLLWithQ( B, QR, t, d, lllCtrl );
        if( ctrl.time )
            Output("Initial LLL time: ",bkz::bkzTimer.Stop()," seconds");
        if( ctrl.progress )
            Output("Initial LLL applied ",lllInfo.numSwaps," swaps");
        numSwaps = lllInfo.numSwaps;
    }
    // The zero columns should be at the end of B
    const Int rank = lllInfo.rank;

    ofstream failedEnumFile, streakSizesFile,
             normsFile, projNormsFile,
             nontrivialCoordsFile;
    if( ctrl.logFailedEnums )
        failedEnumFile.open( ctrl.failedEnumFile.c_str() );
    if( ctrl.logStreakSizes )
        streakSizesFile.open( ctrl.streakSizesFile.c_str() );
    if( ctrl.logNorms )
        normsFile.open( ctrl.normsFile.c_str() );
    if( ctrl.logProjNorms )
        projNormsFile.open( ctrl.projNormsFile.c_str() );
    if( ctrl.logNontrivialCoords )
        nontrivialCoordsFile.open( ctrl.nontrivialCoordsFile.c_str() );

    auto enumCtrl = ctrl.enumCtrl;
    enumCtrl.disablePrecDrop = true;

    Int z=0;
    Int j = ( ctrl.jumpstart ? ctrl.startCol : 0 ) - 1;
    Int numEnums=0, numEnumFailures=0;
    const Int indent = PushIndent(); 
    while( z < rank-1 ) 
    {
        j = Mod(j+1,rank);
        Int bsize = ctrl.blocksize;
        if( ctrl.variableBlocksize )
            bsize = ctrl.blocksizeFunc(j);
        const Int k = Min(j+bsize-1,rank-1);
        const Int h = Min(k+1,rank-1); 
        if( ctrl.checkpoint )
        {
            Write( B, ctrl.checkpointFileBase, ctrl.checkpointFormat, "B" );
            if( j == 0 )
                Write( B, ctrl.tourFileBase, ctrl.checkpointFormat, "B" );
        }
        if( j == 0 )
        {
            if( ctrl.logNorms )
            {
                for( Int j=0; j<n; ++j )
                    normsFile << FrobeniusNorm(B(ALL,IR(j))) << " ";    
                normsFile << endl;
            }
            if( ctrl.logProjNorms )
            {
                for( Int j=0; j<n; ++j )
                    projNormsFile << RealPart(QR(j,j)) << " ";             
                projNormsFile << endl;
            }
        }
        
        Matrix<F> v;
        auto BEnum = B( ALL, IR(j,k+1) );
        auto QREnum = QR( IR(j,k+1), IR(j,k+1) );
        if( ctrl.time )
            bkz::enumTimer.Start();
        if( ctrl.variableEnumType )
            enumCtrl.enumType = ctrl.enumTypeFunc(j);
        const Range<Int> windowInd = IR(j,Min(j+ctrl.multiEnumWindow,k+1));
        auto normUpperBounds = GetRealPartOfDiagonal(QR(windowInd,windowInd));
        Scale( Min(Sqrt(ctrl.lllCtrl.delta),Real(1)), normUpperBounds );
        const auto minPair =
          MultiShortestVectorEnrichment
          ( BEnum, QREnum, normUpperBounds, v, enumCtrl );
        if( ctrl.time )
            Output("Enum/enrich time: ",bkz::enumTimer.Stop()," seconds");
        ++numEnums;

        const Real minProjNorm = minPair.first;
        const Int insertionInd = minPair.second;

        const Real oldProjNorm = RealPart(QREnum(insertionInd,insertionInd));
        const bool keptMin = ( minProjNorm < oldProjNorm );
        if( keptMin )
        {
            if( ctrl.progress )
            {
                Output
                ("Nontrivial enumeration for window of size ",k+1-j,
                 " with j=",j,", z=",z);
                Print( v, "v" );
                Output("insertion index: ",insertionInd);
                Output("oldProjNorm=",oldProjNorm,", minProjNorm=",minProjNorm,", oldProjNorm-minProjNorm=",oldProjNorm-minProjNorm);
            }

            ++numEnumFailures;
            if( ctrl.logFailedEnums )
                failedEnumFile << j << endl;
            if( ctrl.logStreakSizes )
                streakSizesFile << z << endl;
            if( ctrl.logNontrivialCoords )
            {
                for( Int e=0; e<v.Height(); ++e )
                    nontrivialCoordsFile << v(e) << " ";
                nontrivialCoordsFile << endl;
            }
            z = 0;
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
        if( ctrl.time )
            bkz::bkzTimer.Start();
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
            subCtrl.variableBlocksize = false;
            subCtrl.variableEnumType = false;
            subCtrl.recursive = false;
            subCtrl.logFailedEnums = false;
            subCtrl.logStreakSizes = false;
            subCtrl.logNontrivialCoords = false;
            subCtrl.logNorms = false;
            subCtrl.logProjNorms = false;
            subCtrl.checkpoint = false;
            subCtrl.enumCtrl.disablePrecDrop = true;
            subCtrl.enumCtrl.time = false;
            subCtrl.enumCtrl.progress = false;
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
            subLLLCtrl.startCol = ( keptMin ? j : h-1 );
            subLLLCtrl.recursive = false;
            lllInfo = LLLWithQ( BSub, QRSub, tSub, dSub, subLLLCtrl );
            if( lllInfo.numSwaps != 0 )
                changed = true;
            numSwaps += lllInfo.numSwaps;
        }
        if( ctrl.time )
            Output("BKZ time: ",bkz::bkzTimer.Stop()," seconds");
        if( !keptMin )
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
    if( ctrl.logNorms )
        normsFile.close();
    if( ctrl.logProjNorms )
        projNormsFile.close();

    if( ctrl.time )
    {
        Output("Total enumeration time: ",bkz::enumTimer.Total()," seconds");
        Output("Total sub-BKZ time:     ",bkz::bkzTimer.Total()," seconds");
    }

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
    EL_DEBUG_CSE
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

// TODO: Provide a way to display when changing precision without showing
//       all of the backtracks and deep insertions.

namespace bkz {

template<typename RealLower,typename F>
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
    EL_DEBUG_CSE
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

template<typename RealLower,typename F>
bool TryLowerPrecisionMerge
( const Matrix<F>& CL,
  const Matrix<F>& CR,
        Matrix<F>& B,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
        unsigned neededPrec,
        BKZInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( MantissaIsLonger<Base<F>,RealLower>::value &&
        MantissaBits<RealLower>::value >= neededPrec )
    {
        try
        {
            info = LowerPrecisionMerge<RealLower>
              ( CL, CR, B, QR, t, d, ctrl );
            succeeded = true;
        }
        catch( std::exception& e )
        {
            Output("e.what()=",e.what());
        }
    }
    return succeeded;
}

#ifdef EL_HAVE_MPC
template<typename F>
bool TryLowerPrecisionBigFloatMerge
( const Matrix<F>& CL,
  const Matrix<F>& CR,
        Matrix<F>& B,
        Matrix<F>& QR,
        Matrix<F>& t,
        Matrix<Base<F>>& d,
  const BKZCtrl<Base<F>>& ctrl,
        unsigned neededPrec,
        BKZInfo<Base<F>>& info )
{
    bool succeeded = false;
    if( !IsFixedPrecision<Base<F>>::value )
    {
        // Only move down to a lower-precision MPFR type if the jump is
        // substantial. The current value has been naively chosen.
        const mpfr_prec_t minPrecDiff = 32;
        mpfr_prec_t inputPrec = mpfr::Precision();
        if( neededPrec <= inputPrec-minPrecDiff )
        {
            mpfr::SetPrecision( neededPrec );
            try
            {
                info = LowerPrecisionMerge<BigFloat>
                  ( CL, CR, B, QR, t, d, ctrl );
                succeeded = true;
            }
            catch( std::exception& e )
            {
                Output("e.what()=",e.what());
            }
            mpfr::SetPrecision( inputPrec );
        }
    }
    return succeeded;
}
#endif

template<typename F>
BKZInfo<Base<F>>
RecursiveHelper
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  Int numShuffles,
  bool maintainU,
  const BKZCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( maintainU )
        LogicError("Recursive BKZ does not yet support computing U");

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
            Matrix<F> QRL, tL;
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
            Matrix<F> QRR, tR;
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

        bool succeeded = false;
        Int numPrevSwaps = info.numSwaps;
        const bool isInteger = IsInteger( B );
        if( isInteger )
        {
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

            succeeded = TryLowerPrecisionMerge<float>
              ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<double>
                  ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
#ifdef EL_HAVE_QD
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<DoubleDouble>
                  ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<QuadDouble>
                  ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
#elif defined(EL_HAVE_QUAD)
            if( !succeeded )
                succeeded = TryLowerPrecisionMerge<Quad>
                  ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
#endif
#ifdef EL_HAVE_MPC
            if( !succeeded )
                succeeded = TryLowerPrecisionBigFloatMerge
                  ( CL, CR, B, QR, t, d, ctrl, neededPrec, info );
#endif
            if( succeeded )
                info.numSwaps += numPrevSwaps;
        }
        // TODO: Allow for dropping with non-integer vectors?

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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    Matrix<F> R;
    return BKZ( B, R, ctrl );
}

} // namespace El

#endif // ifndef EL_LATTICE_BKZ_HPP
