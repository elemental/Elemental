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
#ifdef EL_HAVE_MPC
        const BigInt n =
          Input("--n","prime lower bound",BigInt(Pow(BigInt(10),BigInt(20))));
        const Int numTrials = Input("--numTrials","num trials",100);

        BigInt p = NextProbablePrime( n );
        for( Int trial=0; trial<numTrials; ++trial ) 
        {
            // Generate a random number between 1 and p
            BigInt a = SampleUniform( BigInt(1), p );
            // Square the random number
            BigInt aSquared = a*a;
            aSquared %= p; 
            // Compute the square-root mod p
            BigInt aSquaredRoot = SqrtModPrime( aSquared, p ); 
            // If a != aSquaredRoot, throw and exception
            if( aSquared != Mod(aSquaredRoot*aSquaredRoot,p) )
            {
                LogicError(aSquaredRoot,"^2 (mod ",p,") != ",a);
            }
        }
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
