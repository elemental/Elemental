/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_MPC
typedef El::BigFloat Real;
#elif defined(EL_HAVE_QUAD)
typedef El::Quad Real;
#else
typedef double Real;
#endif

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const std::string inputBasisFile =
          El::Input("--inputBasisFile","input basis file",
            std::string("../data/number_theory/SVPChallenge40.txt"));
        const bool trans = El::Input("--transpose","transpose input?",true);
        const std::string outputBasisFile =
          El::Input("--outputBasisFile","output basis file",std::string("BKZ"));
        const std::string shortestVecFile =
          El::Input
          ("--shortestVecFile","shortest vector file",std::string("shortest"));
        const Real delta = El::Input("--delta","delta for LLL",Real(0.9999));
        const Real eta =
          El::Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + El::Pow(El::limits::Epsilon<Real>(),Real(0.9)));
        const El::Int varInt =
          El::Input
          ("--variant",
           "0: weak LLL, 1: normal LLL, 2: deep insertion LLL, "
           "3: deep reduction LLL",1);
        const El::Int blocksize =
          El::Input("--blocksize","BKZ blocksize",20);
        const bool variableBsize =
          El::Input("--variableBsize","variable blocksize?",false);
        const bool variableEnumType =
          El::Input("--variableEnumType","variable enum type?",false);
        const El::Int multiEnumWindow =
          El::Input("--multiEnumWindow","window for y-sparse enumeration",15);
        const El::Int phaseLength =
          El::Input("--phaseLength","YSPARSE_ENUM phase length",10);
        const double enqueueProb =
          El::Input("--enqueueProb","enqueue probability?",1.);
        const El::Int progressLevel =
          El::Input("--progressLevel","YSPARSE_ENUM progress level",4);
        const bool presort = El::Input("--presort","presort columns?",false);
        const bool smallestFirst =
          El::Input("--smallestFirst","sort smallest first?",true);
        const bool recursiveLLL =
          El::Input("--recursiveLLL","recursive LLL?",true);
        const bool recursiveBKZ =
          El::Input("--recursiveBKZ","recursive BKZ?",false);
        const El::Int cutoff = El::Input("--cutoff","recursive cutoff",10);
        const bool earlyAbort =
          El::Input("--earlyAbort","early abort BKZ?",false);
        const El::Int numEnumsBeforeAbort =
          El::Input
          ("--numEnumsBeforeAbort","num enums before early aborting",1000);
        const bool subBKZ =
          El::Input
          ("--subBKZ","use BKZ w/ lower blocksize for subproblems?",true);
        const bool subEarlyAbort =
          El::Input("--subEarlyAbort","early abort subproblem?",false);
        const bool jumpstartBKZ =
          El::Input("--jumpstartBKZ","jumpstart BKZ?",false);
        const El::Int startColBKZ =
          El::Input("--startColBKZ","BKZ start column",0);
        const bool timeLLL = El::Input("--timeLLL","time LLL?",false);
        const bool timeBKZ = El::Input("--timeBKZ","time BKZ?",true);
        const bool progressLLL =
          El::Input("--progressLLL","print LLL progress?",false);
        const bool progressBKZ =
          El::Input("--progressBKZ","print BKZ progress?",true);
        const bool print = El::Input("--print","output all matrices?",true);
        const bool logFailedEnums =
          El::Input("--logFailedEnums","log failed enumerations in BKZ?",false);
        const bool logStreakSizes =
          El::Input("--logStreakSizes","log enum streak sizes in BKZ?",false);
        const bool logNontrivialCoords =
          El::Input("--logNontrivialCoords","log nontrivial enum coords?",false);
        const bool logNorms = El::Input("--logNorms","log norms of B?",true);
        const bool logProjNorms =
          El::Input("--logProjNorms","log proj norms of B?",true);
        const bool checkpoint =
          El::Input("--checkpoint","checkpoint each tour?",true);
        const Real targetRatio =
          El::Input("--targetRatio","targeted ratio of GH(L)",Real(1.05));
        const bool timeEnum = El::Input("--timeEnum","time enum?",true);
        const bool innerEnumProgress =
          El::Input("--innerEnumProgress","inner enum progress?",false);
        const bool probEnum =
          El::Input("--probEnum","probabalistic enumeration *after* BKZ?",true);
        const bool fullEnum = El::Input("--fullEnum","SVP via full enum?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec =
          El::Input("--prec","MPFR precision",mpfr_prec_t(1024));
#endif
        El::ProcessInput();
        El::PrintInputReport();

#ifdef EL_HAVE_MPC
        El::mpfr::SetPrecision( prec );
#endif

        El::Matrix<Real> B;
        if( trans )
        {
            El::Matrix<Real> BTrans;
            El::Read( BTrans, inputBasisFile );
            El::Transpose( BTrans, B );
        }
        else
            El::Read( B, inputBasisFile );
        const El::Int m = B.Height();
        const El::Int n = B.Width();
        const Real BOrigOne = El::OneNorm( B );
        El::Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            El::Print( B, "BOrig" );

        auto blocksizeLambda =
          [&]( El::Int j )
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
          [&]( El::Int j )
          {
              if( j <= 3 )
                  return El::YSPARSE_ENUM;
              else
                  return El::FULL_ENUM;
              //return El::FULL_ENUM;
          };
        El::BKZCtrl<Real> ctrl;
        ctrl.blocksize = blocksize;
        ctrl.variableBlocksize = variableBsize;
        ctrl.blocksizeFunc = El::MakeFunction(blocksizeLambda);
        ctrl.variableEnumType = variableEnumType;
        ctrl.enumTypeFunc = El::MakeFunction(enumTypeLambda);
        ctrl.multiEnumWindow = multiEnumWindow;
        ctrl.time = timeBKZ;
        ctrl.progress = progressBKZ;
        ctrl.recursive = recursiveBKZ;
        ctrl.jumpstart = jumpstartBKZ;
        ctrl.startCol = startColBKZ;
        ctrl.enumCtrl.enumType = El::FULL_ENUM;
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
        ctrl.lllCtrl.variant = static_cast<El::LLLVariant>(varInt);
        ctrl.lllCtrl.recursive = recursiveLLL;
        ctrl.lllCtrl.cutoff = cutoff;
        ctrl.lllCtrl.presort = presort;
        ctrl.lllCtrl.smallestFirst = smallestFirst;
        ctrl.lllCtrl.progress = progressLLL;
        ctrl.lllCtrl.time = timeLLL;

        // TODO(poulson): Make this less fragile
        /*
        ctrl.enumCtrl.customMinInfNorms = true;
        ctrl.enumCtrl.customMaxInfNorms = true;
        ctrl.enumCtrl.customMinOneNorms = true;
        ctrl.enumCtrl.customMaxOneNorms = true;
        const El::Int startIndex = Max(n/2-1,0);
        const El::Int numPhases = ((n-startIndex)+phaseLength-1) / phaseLength;
        El::Output("numPhases=",numPhases);
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

        ctrl.enumCtrl.maxOneNorms[0] = 0;
        ctrl.enumCtrl.maxOneNorms[1] = 1;
        ctrl.enumCtrl.maxOneNorms[2] = 1;
        ctrl.enumCtrl.maxOneNorms[3] = 1;
        ctrl.enumCtrl.maxOneNorms[4] = 1;
        ctrl.enumCtrl.maxOneNorms[5] = 2;
        ctrl.enumCtrl.maxOneNorms[6] = 3;
        ctrl.enumCtrl.maxOneNorms[7] = 3;

        ctrl.enumCtrl.maxInfNorms[0] = 1;
        ctrl.enumCtrl.maxInfNorms[1] = 1;
        ctrl.enumCtrl.maxInfNorms[2] = 1;
        ctrl.enumCtrl.maxInfNorms[3] = 1;
        ctrl.enumCtrl.maxInfNorms[4] = 1;
        ctrl.enumCtrl.maxInfNorms[5] = 1;
        ctrl.enumCtrl.maxInfNorms[6] = 2;
        ctrl.enumCtrl.maxInfNorms[7] = 2;
        */

        const double startTime = El::mpi::Time();
        El::Matrix<Real> R;
        auto info = El::BKZ( B, R, ctrl );
        const double runTime = El::mpi::Time() - startTime;
        El::Output
        ("  BKZ(",blocksize,",",delta,",",eta,") took ",runTime," seconds");
        El::Output("    achieved delta:   ",info.delta);
        El::Output("    achieved eta:     ",info.eta);
        El::Output("    num swaps:        ",info.numSwaps);
        El::Output("    num enums:        ",info.numEnums);
        El::Output("    num failed enums: ",info.numEnumFailures);
        El::Output("    log(vol(L)):      ",info.logVol);
        const Real GH = El::LatticeGaussianHeuristic( info.rank, info.logVol );
        const Real challenge = targetRatio*GH;
        El::Output("    GH(L):             ",GH);
        El::Output("    targetRatio*GH(L): ",challenge);
        if( print )
        {
            El::Print( B, "B" );
            El::Print( R, "R" );
        }
        El::Write( B, outputBasisFile, El::ASCII, "BKZ" );
        const Real BOneNorm = El::OneNorm( B );
        El::Output("|| B ||_1 = ",BOneNorm);

        auto b0 = B( El::ALL, El::IR(0) );
        const Real b0Norm = El::FrobeniusNorm( b0 );
        El::Output("|| b_0 ||_2 = ",b0Norm);
        El::Print( b0, "b0" );
        bool succeeded = false;
        if( b0Norm <= challenge )
        {
            El::Output
            ("SVP Challenge solved via BKZ: || b_0 ||_2=",b0Norm,
             " <= targetRatio*GH(L)=",challenge);
            succeeded = true;
            El::Write( b0, shortestVecFile, El::ASCII, "b0" );
        }
        else
            El::Output
            ("SVP Challenge NOT solved via BKZ: || b_0 ||_2=",b0Norm,
             " > targetRatio*GH(L)=",challenge);

        if( !succeeded || fullEnum )
        {
            const El::Int start = 0;
            const El::Int numCols = n;
            const El::Range<El::Int> subInd( start, start+numCols );
            auto BSub = B( El::ALL, subInd );
            auto RSub = R( subInd, subInd );

            const Real target = start == 0 ? challenge : RSub(0,0);

            El::Timer timer;
            El::Matrix<Real> v;
            El::EnumCtrl<Real> enumCtrl;
            enumCtrl.enumType = probEnum ? El::GNR_ENUM : El::FULL_ENUM;
            timer.Start();
            Real result;
            if( fullEnum )
            {
                result =
                  El::ShortestVectorEnumeration
                  ( BSub, RSub, target, v, enumCtrl );
                El::Output("shortest vector result = ",result);
            }
            else
            {
                result =
                  El::ShortVectorEnumeration( BSub, RSub, target, v, enumCtrl );
                El::Output("short vector result = ",result);
            }
            El::Output("Enumeration: ",timer.Stop()," seconds");
            if( result < target )
            {
                El::Print( BSub, "BSub" );
                El::Print( v, "v" );
                El::Matrix<Real> x;
                El::Zeros( x, m, 1 );
                El::Gemv( El::NORMAL, Real(1), BSub, v, Real(0), x );
                El::Print( x, "x" );
                const Real xNorm = El::FrobeniusNorm( x );
                El::Output("|| x ||_2 = ",xNorm);
                El::Output("Claimed || x ||_2 = ",result);
                El::Write( x, shortestVecFile, El::ASCII, "x" );

                El::EnrichLattice( BSub, v );
                El::Print( B, "BNew" );
            }
            else
                El::Output("Enumeration failed after ",timer.Stop()," seconds");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }
    return 0;
}
