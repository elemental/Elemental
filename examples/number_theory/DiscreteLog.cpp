/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        const int minIntBits = Input("--minIntBits","integer bits",384);
        const bool multistage =
          Input("--multistage","Use Pohlig-Hellman?",true);
        const bool assumePrime =
          Input("--assumePrime","assume prime modulus?",false);
        const BigInt a0 = Input("--a0","a0 in Pollard rho",BigInt(0));
        const BigInt b0 = Input("--b0","b0 in Pollard rho",BigInt(0));
        const int numReps = Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = Input("--progress","factor progress?",true);
        const bool time = Input("--time","time Pollard rho steps?",true);
        ProcessInput();
        PrintInputReport(); 

        mpfr::SetMinIntBits( minIntBits );

        dlog::PollardRhoCtrl rhoCtrl;
        rhoCtrl.a0 = a0;
        rhoCtrl.b0 = b0;
        rhoCtrl.multistage = multistage;
        rhoCtrl.assumePrime = assumePrime;
        rhoCtrl.factorCtrl.numReps = numReps;
        rhoCtrl.factorCtrl.progress = progress;
        rhoCtrl.factorCtrl.time = time;
        rhoCtrl.progress = progress;
        rhoCtrl.time = time;

        // p=999959, q=3, r=7
        // We should find an index of 178162
        BigInt q(3), r(7), p(999959);
        Output("Computing discrete log of ",q," w.r.t. ",r," mod ",p);
        BigInt index = dlog::PollardRho( q, r, p, rhoCtrl );
        Output("index: ",index);
        Output("");

        // p=99989, q=107, r=2
        // We should find an index of 87833
        p = 99989;
        q = 107;
        r = 2;
        Output("Computing discrete log of ",q," w.r.t. ",r," mod ",p);
        index = dlog::PollardRho( q, r, p, rhoCtrl );
        Output("index: ",index);
        Output("");

        // Try an impossible problem: 3^x = 4 (mod 13)
        p = 13;
        q = 4;
        r = 3;
        Output("Computing discrete log of ",q," w.r.t. ",r," mod ",p);
        try { index = dlog::PollardRho( q, r, p, rhoCtrl ); }
        catch( std::exception& e ) { Output("exception: ",e.what()); }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
