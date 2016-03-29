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

    try
    {
        typedef unsigned long long TUnsigned;
        const TUnsigned B1=1000000000;
        DynamicSieve<TUnsigned> sieve;
        TUnsigned numPrimes=0;
        for( TUnsigned p=2; p<=B1; p=sieve.NextPrime() )
        {
            //Output(numPrimes,": ",p);
            ++numPrimes; 
        }
        Output("numPrimes=",numPrimes);
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
