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

    try
    {
#ifdef EL_HAVE_MPC
        const El::BigInt n =
          El::Input
          ("--n","prime lower bound",
           El::BigInt(El::Pow(El::BigInt(10),El::BigInt(20))));
        const El::Int numTrials = El::Input("--numTrials","num trials",100);

        El::BigInt p = El::NextProbablePrime( n );
        for( El::Int trial=0; trial<numTrials; ++trial )
        {
            // Generate a random number between 1 and p
            El::BigInt a = El::SampleUniform( El::BigInt(1), p );
            // Square the random number
            El::BigInt aSquared = a*a;
            aSquared %= p;
            // Compute the square-root mod p
            El::BigInt aSquaredRoot = El::SqrtModPrime( aSquared, p );
            // If a != aSquaredRoot, throw and exception
            if( aSquared != El::Mod(aSquaredRoot*aSquaredRoot,p) )
            {
                El::LogicError(aSquaredRoot,"^2 (mod ",p,") != ",a);
            }
        }
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
