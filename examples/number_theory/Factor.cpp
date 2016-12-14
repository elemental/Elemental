/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_MPC
void FactorRho( const El::BigInt& n, const El::factor::PollardRhoCtrl& ctrl )
{
    auto factors = El::factor::PollardRho( n, ctrl );
    El::Output("factors:");
    for( auto factor : factors )
        El::Output("  ",factor);
    El::Output("");
}

template<typename TSieve>
void FactorPM1
( const El::BigInt& n, const El::factor::PollardPMinusOneCtrl<TSieve>& ctrl )
{
    auto factors = El::factor::PollardPMinusOne( n, ctrl );
    El::Output("factors:");
    for( auto factor : factors )
        El::Output("  ",factor);
    El::Output("");
}
#endif

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        typedef unsigned long long TSieve;

        const int minIntBits = El::Input("--minIntBits","integer bits",384);
        const unsigned long numSteps =
          El::Input("--numSteps","x^2k + a in Pollard rho",1u);
        const El::BigInt x0 =
          El::Input("--x0","x0 in Pollard rho",El::BigInt(2));
        const El::Int gcdDelayRho =
          El::Input("--gcdDelayRho","GCD delay in Pollard's rho",100);
        const TSieve smooth1 =
          El::Input
          ("--smooth1","Stage one smoothness bound for (p-1)",1000000ULL);
        const TSieve smooth2 =
          El::Input
          ("--smooth2","Stage two smoothness bound for (p-1)",10000000ULL);
        const El::Int gcdDelayPm1 =
          El::Input("--gcdDelayPm1","GCD delay in Pollard's p-1 (stage 2)",100);
        const bool checkpointPm1 =
          El::Input("--checkpointPm1","checkpoint p-1?",true);
        const El::Int checkpointFreqPm1 =
          El::Input
          ("--checkpointFreqPm1","checkpoint frequency in p-1",1000000);
        const int numReps = El::Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = El::Input("--progress","factor progress?",true);
        const bool time = El::Input("--time","time Pollard rho steps?",true);
        const bool largeRho =
          El::Input("--largeRho","reproduce Pollard and Brent's result?",false);
        const bool heroicPMinusOne =
          El::Input("--heroicPMinusOne","reproduce Zimmerman's result?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::mpfr::SetMinIntBits( minIntBits );

        El::factor::PollardRhoCtrl rhoCtrl;
        rhoCtrl.numSteps = numSteps;
        rhoCtrl.x0 = x0;
        rhoCtrl.gcdDelay = gcdDelayRho;
        rhoCtrl.numReps = numReps;
        rhoCtrl.progress = progress;
        rhoCtrl.time = time;

        El::factor::PollardPMinusOneCtrl<TSieve> pm1Ctrl;
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
        El::BigInt n = El::Pow(El::BigInt(2),unsigned(77)) - 3;
        El::Output("n=2^77-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417)
        n = El::Pow(El::BigInt(2),unsigned(79)) - 3;
        El::Output("n=2^79-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 2^97 - 3
        n = El::Pow(El::BigInt(2),unsigned(97)) - 3;
        El::Output("n=2^97-3=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        // n = 3^100 + 2
        n = El::Pow(El::BigInt(3),unsigned(100)) + 2;
        El::Output("n=3^100+2=",n);
        FactorRho( n, rhoCtrl );
        FactorPM1( n, pm1Ctrl );

        if( largeRho )
        {
            // Try Pollard's rho on Pollard and Brent's famous result of
            //    2^(2^8) + 1 = 1238926361552897 * p_{62},
            // where p_{62} is a 62-digit prime.
            // Note that this should only take about 30 seconds.
            El::mpfr::SetMinIntBits( 1024 );
            El::BigInt z =
              El::Pow(El::BigInt(2),El::Pow(El::BigInt(2),El::BigInt(8))) + 1;
            El::Output("z=2^(2^8)+1=",z);
            FactorRho( z, rhoCtrl );

            if( numSteps != 512u )
            {
                // Try again with numSteps=512 and x0=3
                El::Output("Factoring z=2^(2^8)+1 again with numSteps=512");
                auto rhoCtrlMod = rhoCtrl;
                rhoCtrlMod.numSteps = 512u;
                rhoCtrlMod.x0 = El::BigInt(3);
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
            El::mpfr::SetMinIntBits( 4096 );
            El::BigInt z = El::Pow(El::BigInt(2),El::BigInt(2098)) + 1;
            pm1Ctrl.smooth1 = 5000000000UL;
            pm1Ctrl.smooth2 = 10000000000000UL;
            pm1Ctrl.checkpoint = true;
            pm1Ctrl.checkpointFreq = 1000000;
            El::Output("z=2^2098+1=",z);
            El::DynamicSieve<unsigned long long,unsigned> sieve;
            FactorPM1( z, pm1Ctrl );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }
#endif

    return 0;
}
