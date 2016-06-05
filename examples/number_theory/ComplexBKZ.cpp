/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

typedef double Real;
typedef Complex<Real> F;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","problem size",20);
        const Real radius = Input("--radius","sampling radius",Real(1.e6));
        const Real delta = Input("--delta","delta for LLL",Real(0.9999));
        const Real eta =
          Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9)));
        const Int blocksize = Input("--blocksize","BKZ blocksize",15);
        const bool earlyAbort = Input("--earlyAbort","early abort BKZ?",false);
        const Int numEnumsBeforeAbort =
          Input("--numEnumsBeforeAbort","num enums before early aborting",1000);
        const bool subBKZ =
          Input("--subBKZ","use BKZ w/ lower blocksize for subproblems?",true);
        const bool subEarlyAbort =
          Input("--subEarlyAbort","early abort subproblem?",false);
        const bool timeLLL = Input("--timeLLL","time LLL?",false);
        const bool timeBKZ = Input("--timeBKZ","time BKZ?",true);
        const bool progressLLL =
          Input("--progressLLL","print LLL progress?",false); 
        const bool progressBKZ =
          Input("--progressBKZ","print BKZ progress?",true); 
        const bool print = Input("--print","output all matrices?",true);
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
        ProcessInput();
        PrintInputReport();

        Matrix<F> B; 
        KnapsackTypeBasis( B, n, radius );
        const Int m = B.Height();
        const Real BOrigOne = OneNorm( B ); 
        Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            Print( B, "BOrig" );

        BKZCtrl<Real> ctrl;
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

        const double startTime = mpi::Time();
        Matrix<F> R;
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
        Print( B, "BKZ" );
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
            Print( b0, "b0" );
        }
        else
            Output
            ("SVP Challenge NOT solved via BKZ: || b_0 ||_2=",b0Norm,
             " > targetRatio*GH(L)=",challenge);

        if( !succeeded )
        {
            const Real target = RealPart(R(0,0));

            Timer timer;
            Matrix<F> v;
            EnumCtrl<Real> enumCtrl;
            enumCtrl.enumType = FULL_ENUM;
            timer.Start();
            Real result = 
              ShortestVectorEnumeration( B, R, target, v, enumCtrl );
            Output("Enumeration: ",timer.Stop()," seconds");

            if( result < target )
            {
                Print( v, "v" );
                Matrix<F> x;
                Zeros( x, m, 1 );
                Gemv( NORMAL, F(1), B, v, F(0), x );
                Print( x, "x" );
                const Real xNorm = FrobeniusNorm( x );
                Output("|| x ||_2 = ",xNorm);
                Output("Claimed || x ||_2 = ",result);
                Print( x, "x" );

                EnrichLattice( B, v );
                Print( B, "BFinal" );
            }
            else
                Output("Previous result was optimal");
        }
    }
    catch( std::exception& e ) { ReportException(e); }
    return 0;
}
