/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

#ifdef EL_HAVE_MPC
typedef BigFloat F;
#elif defined(EL_HAVE_QUAD)
typedef Quad F;
#else
typedef double F;
#endif
typedef Base<F> Real;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const string inputBasisFile =
          Input("--inputBasisFile","input basis file",string("SVPChallenge40.txt"));
        const bool trans = Input("--transpose","transpose input?",true);
        const string outputBasisFile = 
          Input("--outputBasisFile","output basis file",string("BKZ"));
        const string shortestVecFile = 
          Input
          ("--shortestVecFile","shortest vector file",string("shortest"));
        const Real delta = Input("--delta","delta for LLL",Real(0.9999));
        const Real eta =
          Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9)));
        const Int varInt = Input("--variant","0: weak LLL, 1: normal LLL, 2: deep insertion LLL, 3: deep reduction LLL",1);
        const Int blocksize = Input("--blocksize","BKZ blocksize",20);
        const bool variableBsize = Input("--variableBsize","variable blocksize?",false);
        const bool variableEnumType = Input("--variableEnumType","variable enum type?",false);
        const Int multiEnumWindow = Input("--multiEnumWindow","window for y-sparse enumeration",15);
        const Int phaseLength =
          Input("--phaseLength","YSPARSE_ENUM phase length",10);
        const double enqueueProb =
          Input("--enqueueProb","enqueue probability?",1.);
        const Int progressLevel =
          Input("--progressLevel","YSPARSE_ENUM progress level",4);
        const bool presort = Input("--presort","presort columns?",false);
        const bool smallestFirst =
          Input("--smallestFirst","sort smallest first?",true);
        const bool recursiveLLL = Input("--recursiveLLL","recursive LLL?",true);
        const bool recursiveBKZ =
          Input("--recursiveBKZ","recursive BKZ?",false);
        const Int cutoff = Input("--cutoff","recursive cutoff",10);
        const bool earlyAbort = Input("--earlyAbort","early abort BKZ?",false);
        const Int numEnumsBeforeAbort =
          Input("--numEnumsBeforeAbort","num enums before early aborting",1000);
        const bool subBKZ =
          Input("--subBKZ","use BKZ w/ lower blocksize for subproblems?",true);
        const bool subEarlyAbort =
          Input("--subEarlyAbort","early abort subproblem?",false);
        const bool jumpstartBKZ =
          Input("--jumpstartBKZ","jumpstart BKZ?",false);
        const Int startColBKZ = Input("--startColBKZ","BKZ start column",0);
        const bool timeLLL = Input("--timeLLL","time LLL?",false);
        const bool timeBKZ = Input("--timeBKZ","time BKZ?",true);
        const bool progressLLL =
          Input("--progressLLL","print LLL progress?",false); 
        const bool progressBKZ =
          Input("--progressBKZ","print BKZ progress?",false); 
        const bool print = Input("--print","output all matrices?",true);
        const bool logFailedEnums =
          Input("--logFailedEnums","log failed enumerations in BKZ?",false);
        const bool logStreakSizes = 
          Input("--logStreakSizes","log enum streak sizes in BKZ?",false);
        const bool logNontrivialCoords =
          Input("--logNontrivialCoords","log nontrivial enum coords?",false);
        const bool logNorms = Input("--logNorms","log norms of B?",true);
        const bool logProjNorms =
          Input("--logProjNorms","log proj norms of B?",true);
        const bool checkpoint =
          Input("--checkpoint","checkpoint each tour?",true);
        const Real targetRatio =
          Input("--targetRatio","targeted ratio of GH(L)",Real(1.05));
        const bool timeEnum = Input("--timeEnum","time enum?",true);
        const bool innerEnumProgress =
          Input("--innerEnumProgress","inner enum progress?",false);
        const bool probEnum =
          Input("--probEnum","probabalistic enumeration *after* BKZ?",true);
        const bool fullEnum = Input("--fullEnum","SVP via full enum?",false);
        const bool enumOnSubset = Input("--enumOnSubset","enum on subset?",false);
        const Int subsetStart = Input("--subsetStart","start of subset",0);
        const Int subsetSize = Input("--subsetSize","num cols in subset",60);
        const bool doubleCycle = Input("--doubleCycle","cycle last vectors?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec =
          Input("--prec","MPFR precision",mpfr_prec_t(1024));
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        Matrix<Real> B; 
        if( trans )
        {
            Matrix<Real> BTrans;
            Read( BTrans, inputBasisFile );
            Transpose( BTrans, B ); 
        }
        else
            Read( B, inputBasisFile );
        const Int m = B.Height();
        const Int n = B.Width();
        const Real BOrigOne = OneNorm( B ); 
        Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            Print( B, "BOrig" );

        auto blocksizeLambda =
          [&]( Int j )
          {
              // With k-sparse
              if( j <= 3 )
                  return 146;
              else if( j <= 10 )
                  return 62;
              else if( j <= 20 )
                  return 60;
              else if( j <= 50 )
                  return 55;
              else
                  return 45;
              // Full enum
              /*
              if( j == 0 )
                  return 80;
              else if( j == 1 )
                  return 75;
              else if( j == 2 )
                  return 70;
              else if( j <= 10 )
                  return 62;
              else if( j <= 20 )
                  return 60;
              else if( j <= 50 )
                  return 55;
              else
                  return 45;
              */
          };
        auto enumTypeLambda = 
          [&]( Int j )
          {
              if( j <= 3 )
                  return YSPARSE_ENUM;
              else
                  return FULL_ENUM;
              //return FULL_ENUM;
          };
        BKZCtrl<Real> ctrl;
        ctrl.blocksize = blocksize;
        ctrl.variableBlocksize = variableBsize;
        ctrl.blocksizeFunc = function<Int(Int)>(blocksizeLambda);
        ctrl.variableEnumType = variableEnumType;
        ctrl.enumTypeFunc = function<EnumType(Int)>(enumTypeLambda);
        ctrl.multiEnumWindow = multiEnumWindow;
        ctrl.time = timeBKZ;
        ctrl.progress = progressBKZ;
        ctrl.recursive = recursiveBKZ;
        ctrl.jumpstart = jumpstartBKZ;
        ctrl.startCol = startColBKZ;
        ctrl.enumCtrl.enumType = FULL_ENUM;
        ctrl.enumCtrl.time = timeEnum;
        ctrl.enumCtrl.innerProgress = innerEnumProgress; 
        ctrl.enumCtrl.phaseLength = phaseLength;
        ctrl.enumCtrl.enqueueProb = enqueueProb;
        ctrl.enumCtrl.progressLevel = progressLevel;
        ctrl.earlyAbort = earlyAbort;
        ctrl.numEnumsBeforeAbort = numEnumsBeforeAbort;
        ctrl.subBKZ = subBKZ;
        ctrl.subEarlyAbort = subEarlyAbort;
        ctrl.logFailedEnums = logFailedEnums;
        ctrl.logStreakSizes = logStreakSizes;
        ctrl.logNontrivialCoords = logNontrivialCoords;
        ctrl.logNorms = logNorms;
        ctrl.logProjNorms = logProjNorms;
        ctrl.checkpoint = checkpoint;
        ctrl.lllCtrl.delta = delta;
        ctrl.lllCtrl.eta = eta;
        ctrl.lllCtrl.variant = static_cast<LLLVariant>(varInt);
        ctrl.lllCtrl.recursive = recursiveLLL;
        ctrl.lllCtrl.cutoff = cutoff;
        ctrl.lllCtrl.presort = presort;
        ctrl.lllCtrl.smallestFirst = smallestFirst;
        ctrl.lllCtrl.progress = progressLLL;
        ctrl.lllCtrl.time = timeLLL;

        ctrl.enumCtrl.customMinInfNorms = true;
        ctrl.enumCtrl.customMaxInfNorms = true;
        ctrl.enumCtrl.customMinOneNorms = true;
        ctrl.enumCtrl.customMaxOneNorms = true;
        const Int startIndex = Max(n/2-1,0);
        const Int numPhases = ((n-startIndex)+phaseLength-1) / phaseLength;
        ctrl.enumCtrl.minInfNorms.resize( numPhases, 0 );
        ctrl.enumCtrl.maxInfNorms.resize( numPhases, 1 );
        ctrl.enumCtrl.minOneNorms.resize( numPhases, 0 );
        ctrl.enumCtrl.maxOneNorms.resize( numPhases, 1 );
        // NOTE: This is tailored to SVP 146 where the ranges are
        // 0: [72,82)
        // 1: [82,92)
        // 2: [92,102)
        // 3: [102,112)
        // 4: [112,122)
        // 5: [122,132)
        // 6: [132,142)
        // 7: [142,146)

        /*
        ctrl.enumCtrl.minOneNorms[0] = 0;
        ctrl.enumCtrl.minOneNorms[1] = 0;
        ctrl.enumCtrl.minOneNorms[2] = 1;
        ctrl.enumCtrl.minOneNorms[3] = 1;
        ctrl.enumCtrl.minOneNorms[4] = 1;
        ctrl.enumCtrl.minOneNorms[5] = 2;
        ctrl.enumCtrl.minOneNorms[6] = 3;
        ctrl.enumCtrl.minOneNorms[7] = 3;
        */

        ctrl.enumCtrl.maxOneNorms[0] = 0;
        ctrl.enumCtrl.maxOneNorms[1] = 1;
        ctrl.enumCtrl.maxOneNorms[2] = 1;
        ctrl.enumCtrl.maxOneNorms[3] = 1;
        ctrl.enumCtrl.maxOneNorms[4] = 1;
        ctrl.enumCtrl.maxOneNorms[5] = 2;
        ctrl.enumCtrl.maxOneNorms[6] = 3;
        ctrl.enumCtrl.maxOneNorms[7] = 3;

        /*
        ctrl.enumCtrl.minInfNorms[0] = 0;
        ctrl.enumCtrl.minInfNorms[1] = 0;
        ctrl.enumCtrl.minInfNorms[2] = 1;
        ctrl.enumCtrl.minInfNorms[3] = 1;
        ctrl.enumCtrl.minInfNorms[4] = 1;
        ctrl.enumCtrl.minInfNorms[5] = 1;
        ctrl.enumCtrl.minInfNorms[6] = 1;
        ctrl.enumCtrl.minInfNorms[7] = 1;
        */

        ctrl.enumCtrl.maxInfNorms[0] = 1;
        ctrl.enumCtrl.maxInfNorms[1] = 1;
        ctrl.enumCtrl.maxInfNorms[2] = 1;
        ctrl.enumCtrl.maxInfNorms[3] = 1;
        ctrl.enumCtrl.maxInfNorms[4] = 1;
        ctrl.enumCtrl.maxInfNorms[5] = 1;
        ctrl.enumCtrl.maxInfNorms[6] = 2;
        ctrl.enumCtrl.maxInfNorms[7] = 2;

        const double startTime = mpi::Time();
        Matrix<Real> R;
        auto info = BKZ( B, R, ctrl );
        const double runTime = mpi::Time() - startTime;
        Output
        ("  BKZ(",blocksize,",",delta,",",eta,") took ",runTime," seconds"); 
        Output("    achieved delta:   ",info.delta);
        Output("    achieved eta:     ",info.eta);
        Output("    num swaps:        ",info.numSwaps);
        Output("    num enums:        ",info.numEnums);
        Output("    num failed enums: ",info.numEnumFailures);
        Output("    log(vol(L)):      ",info.logVol);
        const Real GH = LatticeGaussianHeuristic( info.rank, info.logVol );
        const Real challenge = targetRatio*GH;
        Output("    GH(L):             ",GH);
        Output("    targetRatio*GH(L): ",challenge);
        if( print )
        {
            Print( B, "B" ); 
            Print( R, "R" );
        }
        Write( B, outputBasisFile, ASCII, "BKZ" );
        const Real BOneNorm = OneNorm( B );
        Output("|| B ||_1 = ",BOneNorm);

        auto b0 = B( ALL, IR(0) );
        const Real b0Norm = FrobeniusNorm( b0 );
        Output("|| b_0 ||_2 = ",b0Norm);
        if( print )
            Print( b0, "b0" );
        bool succeeded = false;
        if( b0Norm <= challenge )
        {
            Output
            ("SVP Challenge solved via BKZ: || b_0 ||_2=",b0Norm,
             " <= targetRatio*GH(L)=",challenge);
            succeeded = true;
            Write( b0, shortestVecFile, ASCII, "b0" );
        }
        else
            Output
            ("SVP Challenge NOT solved via BKZ: || b_0 ||_2=",b0Norm,
             " > targetRatio*GH(L)=",challenge);

        if( !succeeded || fullEnum || (enumOnSubset && subsetStart != 0) )
        {
            const Int start = ( enumOnSubset ? subsetStart : 0 ); 
            const Int numCols = ( enumOnSubset ? subsetSize : n );
            const Range<Int> subInd( start, start+numCols );
            auto BSub = B( ALL, subInd );
            auto RSub = R( subInd, subInd );

            const Real target = ( start == 0 ? challenge : RSub.Get(0,0) );

            Timer timer;
            if( enumOnSubset && doubleCycle && subsetSize >= 2 )
            {
                Matrix<double> v;
                EnumCtrl<double> enumCtrl;
                enumCtrl.enumType = ( probEnum ? GNR_ENUM : FULL_ENUM );
                enumCtrl.numTrials = 1;

                Matrix<double> BSubSwap;
                Zeros( BSubSwap, m, subsetSize );
                auto BL = B( ALL, IR(start,start+subsetSize-2) );
                auto BSubSwapL = BSubSwap( ALL, IR(0,subsetSize-2) );
                Copy( BL, BSubSwapL );
                for( Int j=start+subsetSize-2; j<n-1; ++j )
                {
                    auto bj = B( ALL, IR(j) ); 
                    auto bSubSwapj = BSubSwap( ALL, IR(subsetSize-2) );
                    Copy( bj, bSubSwapj );
                    for( Int k=j+1; k<n-1; ++k )
                    {
                        auto bk = B( ALL, IR(k) );
                        auto bSubSwapk = BSubSwap( ALL, IR(subsetSize-1) ); 
                        Copy( bk, bSubSwapk );
                        Matrix<double> RSubSwap( BSubSwap );
                        Output("Cycling with j=",j,", k=",k);
                        qr::ExplicitTriang( RSubSwap );
                        timer.Start();
                        Real result =
                          ShortestVectorEnumeration
                          ( BSubSwap, RSubSwap, double(target), v, enumCtrl );
                        Output("Enumeration: ",timer.Stop()," seconds");
                        if( result < RSubSwap.Get(0,0)-double(0.001) )
                        {
                            Print( BSubSwap, "BSubSwap" );
                            Print( v, "v" );
                            Matrix<double> x;
                            Zeros( x, m, 1 );
                            Gemv( NORMAL, 1., BSubSwap, v, 0., x );
                            Print( x, "x" );
                            const double xNorm = FrobeniusNorm( x );
                            Output("|| x ||_2 = ",xNorm);
                            Output("Claimed || x ||_2 = ",result);
                            Write( x, shortestVecFile, ASCII, "x" );
                        }
                    }
                }
            }
            else
            {
                Matrix<F> v;
                EnumCtrl<Real> enumCtrl;
                enumCtrl.enumType = ( probEnum ? GNR_ENUM : FULL_ENUM );
                timer.Start();
                Real result;
                if( fullEnum )
                  result = 
                    ShortestVectorEnumeration( BSub, RSub, target, v, enumCtrl );
                else
                  result = 
                    ShortVectorEnumeration( BSub, RSub, target, v, enumCtrl );
                Output("Enumeration: ",timer.Stop()," seconds");
                if( result < target )
                {
                    Print( BSub, "BSub" );
                    Print( v, "v" );
                    Matrix<Real> x;
                    Zeros( x, m, 1 );
                    Gemv( NORMAL, Real(1), BSub, v, Real(0), x );
                    Print( x, "x" );
                    const Real xNorm = FrobeniusNorm( x );
                    Output("|| x ||_2 = ",xNorm);
                    Output("Claimed || x ||_2 = ",result);
                    Write( x, shortestVecFile, ASCII, "x" );

                    EnrichLattice( BSub, v );
                    Print( B, "BNew" );
                }
                else
                    Output("Enumeration failed after ",timer.Stop()," seconds");
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }
    return 0;
}
