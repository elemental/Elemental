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
        const string outputBasisFile = 
          Input("--outputBasisFile","output basis file",string("BKZ"));
        const string shortestVecFile = 
          Input
          ("--shortestVecFile","shortest vector file",string("shortest"));
        const Real delta = Input("--delta","delta for LLL",Real(0.99));
        const Real eta =
          Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9)));
        const Int varInt = Input("--variant","0: weak LLL, 1: normal LLL, 2: deep insertion LLL, 3: deep reduction LLL",1);
        const Int blocksize = Input("--blocksize","BKZ blocksize",20);
        const bool presort = Input("--presort","presort columns?",false);
        const bool smallestFirst =
          Input("--smallestFirst","sort smallest first?",true);
        const bool recursiveLLL = Input("--recursiveLLL","recursive LLL?",true);
        const bool recursiveBKZ =
          Input("--recursiveBKZ","recursive BKZ?",false);
        const Int cutoff = Input("--cutoff","recursive cutoff",10);
        const bool earlyAbort = Input("--earlyAbort","early abort BKZ?",false);
        const bool subBKZ =
          Input("--subBKZ","use BKZ w/ lower blocksize for subproblems?",true);
        const bool subEarlyAbort =
          Input("--subEarlyAbort","early abort subproblem?",false);
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
        const Real targetRatio =
          Input("--targetRatio","targeted ratio of GH(L)",Real(1.05));
        const bool probBKZEnum =
          Input("--probBKZEnum","probabalistic enumeration *within* BKZ?",false);
        const bool probEnum =
          Input("--probEnum","probabalistic enumeration *after* BKZ?",true);
        const bool fullEnum = Input("--fullEnum","SVP via full enum?",false);
        const bool trans = Input("--transpose","transpose input?",true);
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
        const Real BOrigOne = OneNorm( B ); 
        Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            Print( B, "BOrig" );

        BKZCtrl<Real> ctrl;
        ctrl.blocksize = blocksize;
        ctrl.time = timeBKZ;
        ctrl.progress = progressBKZ;
        ctrl.recursive = recursiveBKZ;
        ctrl.enumCtrl.probabalistic = probBKZEnum;
        ctrl.earlyAbort = earlyAbort;
        ctrl.subBKZ = subBKZ;
        ctrl.subEarlyAbort = subEarlyAbort;
        ctrl.logFailedEnums = logFailedEnums;
        ctrl.logStreakSizes = logStreakSizes;
        ctrl.logNontrivialCoords = logNontrivialCoords;
        ctrl.lllCtrl.delta = delta;
        ctrl.lllCtrl.eta = eta;
        ctrl.lllCtrl.variant = static_cast<LLLVariant>(varInt);
        ctrl.lllCtrl.recursive = recursiveLLL;
        ctrl.lllCtrl.cutoff = cutoff;
        ctrl.lllCtrl.presort = presort;
        ctrl.lllCtrl.smallestFirst = smallestFirst;
        ctrl.lllCtrl.progress = progressLLL;
        ctrl.lllCtrl.time = timeLLL;

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

        if( !succeeded || fullEnum )
        {
            // TODO: Provide support for enumerating over the first, say,
            //       2/3 of the vectors
            Matrix<F> v;
            Timer timer;
            timer.Start();
            Real result;
            EnumCtrl<Real> enumCtrl;
            enumCtrl.probabalistic = probEnum;
            if( fullEnum )
              result = 
                ShortestVectorEnumeration( B, R, challenge, v, enumCtrl );
            else
              result = 
                ShortVectorEnumeration( B, R, challenge, v, enumCtrl );
            if( result < challenge )
            {
                Output("Enumeration: ",timer.Stop()," seconds");
                Print( B, "B" );
                Print( v, "v" );
                const Int m = B.Height();
                Matrix<Real> x;
                Zeros( x, m, 1 );
                Gemv( NORMAL, Real(1), B, v, Real(0), x );
                Print( x, "x" );
                const Real xNorm = FrobeniusNorm( x );
                Output("|| x ||_2 = ",xNorm);
                Output("Claimed || x ||_2 = ",result);
                Write( x, shortestVecFile, ASCII, "x" );
            }
            else
                Output("Enumeration failed after ",timer.Stop()," seconds");
        }
    }
    catch( std::exception& e ) { ReportException(e); }
    return 0;
}
