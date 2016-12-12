/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int n = El::Input("--n","problem size",20);
        const Real radius = El::Input("--radius","sampling radius",Real(1.e6));
        const Real delta = El::Input("--delta","delta for LLL",Real(0.9999));
        const Real eta =
          El::Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + El::Pow(El::limits::Epsilon<Real>(),Real(0.9)));
        const El::Int blocksize = El::Input("--blocksize","BKZ blocksize",15);
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
        const bool timeLLL = El::Input("--timeLLL","time LLL?",false);
        const bool timeBKZ = El::Input("--timeBKZ","time BKZ?",true);
        const bool progressLLL =
          El::Input("--progressLLL","print LLL progress?",false);
        const bool progressBKZ =
          El::Input("--progressBKZ","print BKZ progress?",true);
        const bool print = El::Input("--print","output all matrices?",true);
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
        El::ProcessInput();
        El::PrintInputReport();

        El::Matrix<Scalar> B;
        El::KnapsackTypeBasis( B, n, radius );
        const Real BOrigOne = El::OneNorm( B );
        El::Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            El::Print( B, "BOrig" );

        El::BKZCtrl<Real> ctrl;
        ctrl.blocksize = blocksize;
        ctrl.time = timeBKZ;
        ctrl.progress = progressBKZ;
        ctrl.enumCtrl.time = timeEnum;
        ctrl.enumCtrl.innerProgress = innerEnumProgress;
        ctrl.earlyAbort = earlyAbort;
        ctrl.numEnumsBeforeAbort = numEnumsBeforeAbort;
        ctrl.subBKZ = subBKZ;
        ctrl.subEarlyAbort = subEarlyAbort;
        ctrl.logNorms = logNorms;
        ctrl.logProjNorms = logProjNorms;
        ctrl.checkpoint = checkpoint;
        ctrl.lllCtrl.delta = delta;
        ctrl.lllCtrl.eta = eta;
        ctrl.lllCtrl.progress = progressLLL;
        ctrl.lllCtrl.time = timeLLL;

        const double startTime = El::mpi::Time();
        El::Matrix<Scalar> R;
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
        const Real BOneNorm = El::OneNorm( B );
        El::Output("|| B ||_1 = ",BOneNorm);

        auto b0 = B( El::ALL, El::IR(0) );
        const Real b0Norm = El::FrobeniusNorm( b0 );
        El::Output("|| b_0 ||_2 = ",b0Norm);
        El::Print( b0, "b0" );
    }
    catch( std::exception& e ) { El::ReportException(e); }
    return 0;
}
