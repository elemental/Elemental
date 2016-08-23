/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

#ifdef EL_HAVE_MPC
void FactorRho( const BigInt& n, const factor::PollardRhoCtrl& ctrl )
{
    auto factors = factor::PollardRho( n, ctrl );
    Output("factors:");
    for( auto factor : factors )
        Output("  ",factor); 
    Output("");
}

template<typename TSieve>
void FactorPM1
( const BigInt& n, const factor::PollardPMinusOneCtrl<TSieve>& ctrl )
{
    auto factors = factor::PollardPMinusOne( n, ctrl );
    Output("factors:");
    for( auto factor : factors )
        Output("  ",factor); 
    Output("");
}
#endif

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        typedef unsigned long long TSieve;

        const int minIntBits = Input("--minIntBits","integer bits",384);
        const unsigned long numSteps =
          Input("--numSteps","x^2k + a in Pollard rho",1u);
        const BigInt x0 = Input("--x0","x0 in Pollard rho",BigInt(2));
        const Int gcdDelayRho =
          Input("--gcdDelayRho","GCD delay in Pollard's rho",100);
        const TSieve smooth1 =
          Input("--smooth1","Stage one smoothness bound for (p-1)",1000000ULL);
        const TSieve smooth2 =
          Input("--smooth2","Stage two smoothness bound for (p-1)",10000000ULL);
        const Int gcdDelayPm1 = 
          Input("--gcdDelayPm1","GCD delay in Pollard's p-1 (stage 2)",100);
        const bool checkpointPm1 =
          Input("--checkpointPm1","checkpoint p-1?",true);
        const Int checkpointFreqPm1 =
          Input("--checkpointFreqPm1","checkpoint frequency in p-1",1000000);
        const int numReps = Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = Input("--progress","factor progress?",true);
        const bool time = Input("--time","time Pollard rho steps?",true);
        const bool largeRho =
          Input("--largeRho","reproduce Pollard and Brent's result?",false);
        const bool heroicPMinusOne =
          Input("--heroicPMinusOne","reproduce Zimmerman's result?",false);
        ProcessInput();
        PrintInputReport(); 

        mpfr::SetMinIntBits( minIntBits );

        factor::PollardRhoCtrl rhoCtrl;
        rhoCtrl.numSteps = numSteps;
        rhoCtrl.x0 = x0;
        rhoCtrl.gcdDelay = gcdDelayRho;
        rhoCtrl.numReps = numReps;
        rhoCtrl.progress = progress;
        rhoCtrl.time = time;

        factor::PollardPMinusOneCtrl<TSieve> pm1Ctrl;
        pm1Ctrl.smooth1 = smooth1;
        pm1Ctrl.smooth2 = smooth2;
        pm1Ctrl.gcdDelay2 = gcdDelayPm1;
        pm1Ctrl.numReps = numReps;
        pm1Ctrl.progress = progress;
        pm1Ctrl.time = time;
        pm1Ctrl.checkpoint = checkpointPm1;
        pm1Ctrl.checkpointFreq = checkpointFreqPm1;

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) 
        BigInt n = Pow(BigInt(2),unsigned(77)) - 3;
        Output("n=2^77-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417)
        n = Pow(BigInt(2),unsigned(79)) - 3;
        Output("n=2^79-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 2^97 - 3
        n = Pow(BigInt(2),unsigned(97)) - 3;
        Output("n=2^97-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 3^100 + 2 
        n = Pow(BigInt(3),unsigned(100)) + 2;
        Output("n=3^100+2=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        if( largeRho )
        {
            // Try Pollard's rho on Pollard and Brent's famous result of
            //    2^(2^8) + 1 = 1238926361552897 * p_{62},
            // where p_{62} is a 62-digit prime.
            // Note that this should only take about 30 seconds.
            mpfr::SetMinIntBits( 1024 );
            BigInt z = Pow(BigInt(2),Pow(BigInt(2),BigInt(8))) + 1; 
            Output("z=2^(2^8)+1=",z); 
            FactorRho( z, rhoCtrl );

            if( numSteps != 512u )
            {
                // Try again with numSteps=512 and x0=3
                Output("Factoring z=2^(2^8)+1 again with numSteps=512");
                auto rhoCtrlMod = rhoCtrl;
                rhoCtrlMod.numSteps = 512u;
                rhoCtrlMod.x0 = BigInt(3);
                FactorRho( z, rhoCtrlMod );
            }
        }

        if( heroicPMinusOne )
        {
            // Try Pollard's p-1 on 2^2098+1, which P. Zimmerman showed to have
            // a factor of
            // p=372098406910139347411473978297737029649599583843164650153,
            // where p-1=23*32*1049*1627*139999*1284223*7475317*341342347*
            //           2456044907*9909876848747
            mpfr::SetMinIntBits( 4096 );
            BigInt z = Pow(BigInt(2),BigInt(2098)) + 1;
            pm1Ctrl.smooth1 = 5*Pow(10UL,9UL);
            pm1Ctrl.smooth2 = Pow(10UL,13UL);
            pm1Ctrl.checkpoint = true;
            pm1Ctrl.checkpointFreq = 1000000;
            Output("z=2^2098+1=",z);
            DynamicSieve<unsigned long long,unsigned> sieve;
            FactorPM1( z, pm1Ctrl );
        }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
