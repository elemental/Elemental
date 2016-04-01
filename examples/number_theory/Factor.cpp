/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        typedef unsigned long long TSieve;
        typedef unsigned TSieveSmall;

        const int minIntBits = Input("--minIntBits","integer bits",384);
        const unsigned long numSteps =
          Input("--numSteps","x^2k + a in Pollard rho",1u);
        const BigInt x0 = Input("--x0","x0 in Pollard rho",BigInt(2));
        const Int gcdDelay =
          Input("--gcdDelay","GCD delay in Pollard's rho",100);
        const int numReps = Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = Input("--progress","factor progress?",true);
        const bool time = Input("--time","time Pollard rho steps?",true);
        const bool largeRho =
          Input("--largeRho","reproduce Pollard and Brent's result?",false);
        const bool heroicPMinusOne =
          Input("--heroicPMinusOne","reproduce Zimmerman's result?",false);
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        factor::PollardRhoCtrl rhoCtrl;
        rhoCtrl.numSteps = numSteps;
        rhoCtrl.x0 = x0;
        rhoCtrl.gcdDelay = gcdDelay;
        rhoCtrl.numReps = numReps;
        rhoCtrl.progress = progress;
        rhoCtrl.time = time;

        factor::PollardPMinusOneCtrl<TSieve> pm1Ctrl;
        pm1Ctrl.numReps = numReps;
        pm1Ctrl.progress = progress;
        pm1Ctrl.time = time;

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) 
        BigInt n = Pow(BigInt(2),unsigned(77)) - 3;
        Output("n=2^77-3=",n);
        auto factors = factor::PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417)
        n = Pow(BigInt(2),unsigned(79)) - 3;
        Output("n=2^79-3=",n);
        factors = factor::PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 2^97 - 3
        n = Pow(BigInt(2),unsigned(97)) - 3;
        Output("n=2^97-3=",n);
        factors = factor::PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 3^100 + 2 
        n = Pow(BigInt(3),unsigned(100)) + 2;
        Output("n=3^100+2=",n);
        factors = factor::PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        if( largeRho )
        {
            // Try Pollard's rho on Pollard and Brent's famous result of
            //    2^(2^8) + 1 = 1238926361552897 * p_{62},
            // where p_{62} is a 62-digit prime.
            // Note that this should only take about 30 seconds.
            mpc::SetMinIntBits( 1024 );
            BigInt z = Pow(BigInt(2),Pow(BigInt(2),BigInt(8))) + 1; 
            Output("z=2^(2^8)+1=",z); 
            factors = factor::PollardRho( z, rhoCtrl );
            Output("factors:");
            for( auto factor : factors )
                Output("  ",factor); 
            Output("");

            if( numSteps != 512u )
            {
                // Try again with numSteps=512 and x0=3
                Output("Factoring z=2^(2^8)+1 again with numSteps=512");
                auto rhoCtrlMod = rhoCtrl;
                rhoCtrlMod.numSteps = 512u;
                rhoCtrlMod.x0 = BigInt(3);
                factors = factor::PollardRho( z, rhoCtrlMod );
                Output("factors:");
                for( auto factor : factors )
                    Output("  ",factor); 
                Output("");
            }
        }

        if( heroicPMinusOne )
        {
            // Try Pollard's p-1 on 2^2098+1, which P. Zimmerman showed to have
            // a factor of
            // p=372098406910139347411473978297737029649599583843164650153,
            // where p-1=23*32*1049*1627*139999*1284223*7475317*341342347*
            //           2456044907*9909876848747
            mpc::SetMinIntBits( 8192 );
            BigInt z = Pow(BigInt(2),BigInt(2098)) + 1;
            pm1Ctrl.smooth1 = 5*Pow(10UL,9UL);
            pm1Ctrl.smooth2 = Pow(10UL,13UL);
            Output("z=2^2098+1=",z);
            factors = factor::PollardPMinusOne( z, pm1Ctrl );
            Output("factors:");
            for( auto factor : factors )
                Output("  ",factor); 
        }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
