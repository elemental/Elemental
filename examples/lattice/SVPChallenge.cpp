/*
   Copyright (c) 2009-2015, Jack Poulson
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
        const string filename =
          Input("--filename","input file",string("SVPChallenge40.txt"));
        const Real delta = Input("--delta","delta for LLL",Real(0.99));
        const Real eta =
          Input
          ("--eta","eta for LLL",
           Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9)));
        const Int varInt = Input("--variant","0: weak LLL, 1: normal LLL, 2: deep insertion LLL, 3: deep reduction LLL",2);
        const Int blocksize = Input("--blocksize","BKZ blocksize",20);
        const bool presort = Input("--presort","presort columns?",false);
        const bool smallestFirst =
          Input("--smallestFirst","sort smallest first?",true);
        const bool recursive = Input("--recursive","try recursive LLL?",true);
        const Int cutoff = Input("--cutoff","recursive cutoff",10);
        const bool progress = Input("--progress","print progress?",false); 
        const bool time = Input("--time","time LLL?",false);
        const bool print = Input("--print","output all matrices?",true);
        const Real targetRatio =
          Input("--targetRatio","targeted ratio of GH(L)",Real(1.05));
        const bool probEnum =
          Input("--probEnum","probabalistic enumeration?",true);
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
            Read( BTrans, filename );
            Transpose( BTrans, B ); 
        }
        else
            Read( B, filename );
        const Real BOrigOne = OneNorm( B ); 
        Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            Print( B, "BOrig" );

        LLLCtrl<Real> lllCtrl;
        lllCtrl.delta = delta;
        lllCtrl.eta = eta;
        lllCtrl.variant = static_cast<LLLVariant>(varInt);
        lllCtrl.presort = presort;
        lllCtrl.smallestFirst = smallestFirst;
        lllCtrl.progress = progress;
        lllCtrl.time = time;

        BKZCtrl<Real> ctrl;
        ctrl.blocksize = blocksize;
        ctrl.lllCtrl = lllCtrl;

        const double startTime = mpi::Time();
        BKZInfo<Real> info;
        Matrix<Real> R;
        if( recursive )
            info = RecursiveBKZ( B, R, cutoff, ctrl );
        else
            info = BKZ( B, R, ctrl );
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
        }
        else
            Output
            ("SVP Challenge NOT solved via BKZ: || b_0 ||_2=",b0Norm,
             " > targetRatio*GH(L)=",challenge);

        if( !succeeded )
        {
            Matrix<F> v;
            Timer timer;
            timer.Start();
            Real result;
            if( fullEnum )
              result = 
                ShortestVectorEnumeration( B, R, challenge, v, probEnum );
            else
              result = 
                ShortVectorEnumeration( B, R, challenge, v, probEnum );
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
            }
            else
                Output("Enumeration failed after ",timer.Stop()," seconds");
        }
    }
    catch( std::exception& e ) { ReportException(e); }
    return 0;
}
