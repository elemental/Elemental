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
        const Int varInt = Input("--variant","0: weak, 1: normal, 2: deep insertion, 3: deep reduction",2);
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
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec =
          Input("--prec","MPFR precision",mpfr_prec_t(1024));
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        Matrix<Real> BTrans, B;
        Read( BTrans, filename );
        Transpose( BTrans, B ); 
        const Real BOrigOne = OneNorm( B ); 
        Output("|| B_orig ||_1 = ",BOrigOne);
        if( print )
            Print( B, "BOrig" );

        // For now, we simply use LLL to attempt to solve the SVP.
        // Obviously, this will not be competitive before making use
        // of Fincke-Pohst
        LLLCtrl<Real> ctrl;
        ctrl.delta = delta;
        ctrl.eta = eta;
        ctrl.variant = static_cast<LLLVariant>(varInt);
        ctrl.presort = presort;
        ctrl.smallestFirst = smallestFirst;
        ctrl.progress = progress;
        ctrl.time = time;
        const double startTime = mpi::Time();
        LLLInfo<Real> info;
        Matrix<Real> R;
        if( recursive )
            info = RecursiveLLL( B, R, cutoff, ctrl );
        else
            info = LLL( B, R, ctrl );
        const double runTime = mpi::Time() - startTime;
        Output("  LLL(",delta,",",eta,") took ",runTime," seconds"); 
        Output("    achieved delta: ",info.delta);
        Output("    achieved eta:   ",info.eta);
        Output("    num swaps:      ",info.numSwaps);
        Output("    log(vol(L)):    ",info.logVol);
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
            ("SVP Challenge solved via LLL: || b_0 ||_2=",b0Norm,
             " <= targetRatio*GH(L)=",challenge);
            succeeded = true;
        }
        else
            Output
            ("SVP Challenge NOT solved via LLL: || b_0 ||_2=",b0Norm,
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
