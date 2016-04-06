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
    typedef unsigned long long TSieve;
    typedef unsigned TSieveSmall;

    Environment env( argc, argv );
    const TSieve B1 = Input("--B1","smoothness limit",TSieve(1000000000UL));
    const bool print = Input("--print","print primes?",false);
    ProcessInput();
    PrintInputReport();

    try
    {
        Timer timer;

        // Generate all of the primes below the given bound at once
        {
            timer.Start();

            DynamicSieve<TSieve,TSieveSmall> sieve;
            sieve.Generate( B1 );

            Output
            ("Generated all primes up to ",sieve.oddPrimes.back()," in ",
             timer.Stop()," seconds");

            Output("numPrimes (extra): ",sieve.oddPrimes.size()+1);
            while( sieve.oddPrimes.back() > B1 )
                sieve.oddPrimes.pop_back();
            Output("numPrimes=",sieve.oddPrimes.size()+1);

            if( print )
            {
                TSieve p=2;
                for( TSieve i=0; i<sieve.oddPrimes.size(); ++i )
                {
                    Output(i,": ",p);
                    p = sieve.oddPrimes[i];
                }
                Output(sieve.oddPrimes.size(),": ",p);
            }
        }

        // Count the number of primes below the given bound by sequentially
        // generating each
        {
            timer.Start();

            DynamicSieve<TSieve,TSieveSmall> sieve;
            TSieve numPrimes=0;
            for( TSieve p=2; p<=B1; p=sieve.NextPrime() )
            {
                if( print )
                    Output(numPrimes,": ",p);
                ++numPrimes; 
            }

            Output
            ("Iterated over primes below ",B1," in ",timer.Stop()," seconds");
            Output("numPrimes=",numPrimes);
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
