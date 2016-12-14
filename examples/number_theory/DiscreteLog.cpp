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

#ifdef EL_HAVE_MPC
    try
    {
        const int minIntBits = El::Input("--minIntBits","integer bits",384);
        const bool multistage =
          El::Input("--multistage","Use Pohlig-Hellman?",true);
        const bool assumePrime =
          El::Input("--assumePrime","assume prime modulus?",false);
        const El::BigInt a0 =
          El::Input("--a0","a0 in Pollard rho",El::BigInt(0));
        const El::BigInt b0 =
          El::Input("--b0","b0 in Pollard rho",El::BigInt(0));
        const int numReps = El::Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = El::Input("--progress","factor progress?",true);
        const bool time = El::Input("--time","time Pollard rho steps?",true);
        El::ProcessInput();
        El::PrintInputReport();

        El::mpfr::SetMinIntBits( minIntBits );

        El::dlog::PollardRhoCtrl rhoCtrl;
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
        El::BigInt q(3), r(7), p(999959);
        El::Output("Computing discrete log of ",q," w.r.t. ",r," mod ",p);
        El::BigInt index = El::dlog::PollardRho( q, r, p, rhoCtrl );
        El::Output("index: ",index);
        El::Output("");

        // p=99989, q=107, r=2
        // We should find an index of 87833
        p = 99989;
        q = 107;
        r = 2;
        El::Output("Computing discrete log of ",q," w.r.t. ",r," mod ",p);
        index = El::dlog::PollardRho( q, r, p, rhoCtrl );
        El::Output("index: ",index);
        El::Output("");

        // Try an impossible problem: 3^x = 4 (mod 13)
        p = 13;
        q = 4;
        r = 3;
        El::Output
        ("Trying infeasible discrete log of ",q," w.r.t. ",r," mod ",p);
        try { index = El::dlog::PollardRho( q, r, p, rhoCtrl ); }
        catch( std::exception& e ) { El::Output("exception: ",e.what()); }
    }
    catch( std::exception& e ) { El::ReportException(e); }
#endif

    return 0;
}
