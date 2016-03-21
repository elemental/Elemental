/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

#ifdef EL_HAVE_MPC
// x_{i+1} := x_i^2 + a
BigInt PollardRhoBounded
( const BigInt& n, Int a=1, Int gcdDelay=100, Int maxIts=10000 )
{
    BigInt xi=2, x2i=2, Qi=1;
    BigInt xiSave=2, x2iSave=2;
    for( Int i=1; i<=maxIts; ++i )
    {
        // Advance xi once
        xi = Mod( xi*xi + a, n );

        // Advance x2i twice
        x2i = Mod( x2i*x2i + a, n );
        x2i = Mod( x2i*x2i + a, n );

        // Advance Qi
        Qi = Mod( Qi*(x2i-xi), n );

        if( i % gcdDelay == 0 )
        {
            BigInt d = GCD( Qi, n );
            if( d > BigInt(1) )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( d == BigInt(n) )
                {
                    if( gcdDelay == 1 )
                    {
                        Output("(x_n) converged before (x_n mod p) at i=",i);
                        return d;
                    }
                    else
                    {
                        Output("Backtracking and setting gcdDelay=1 at i=",i);
                        i = Max( i-(gcdDelay+1), 0 );
                        gcdDelay = 1;
                        xi = xiSave;
                        x2i = x2iSave;
                    }
                }
                else
                {
                    Output("Found factor d=",d," at i=",i); 
                    return d;
                }
            }

            // NOTE: This was not suggested by Pollard's original paper
            Qi = 1;
            xiSave = xi;
            x2iSave = x2i;
        }
        if( i == maxIts )
        {
            Output("Failed to find a factor in time");
            Output("Final xi=",xi,", x2i=",x2i,", and Qi=",Qi);
        }
    }
    return n;
}

// x_{i+1} := x_i^2 + a
BigInt PollardRho( const BigInt& n, Int a=1, Int gcdDelay=100 )
{
    BigInt xi=2, x2i=2, Qi=1;
    BigInt xiSave=xi, x2iSave=x2i;
    Int k=1, i=1; // it is okay for i to overflow since it is just for printing
    while( true )
    {
        // Advance xi once
        xi = Mod( xi*xi + a, n );

        // Advance x2i twice
        x2i = Mod( x2i*x2i + a, n );
        x2i = Mod( x2i*x2i + a, n );

        // Advance Qi
        Qi = Mod( Qi*(x2i-xi), n );

        if( k >= gcdDelay )
        {
            BigInt d = GCD( Qi, n );
            if( d > BigInt(1) )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( d == BigInt(n) )
                {
                    if( gcdDelay == 1 )
                    {
                        Output("(x_n) converged before (x_n mod p) at i=",i);
                        return d;
                    }
                    else
                    {
                        Output("Backtracking and setting gcdDelay=1 at i=",i);
                        i = Max( i-(gcdDelay+1), 0 );
                        gcdDelay = 1;
                        xi = xiSave;
                        x2i = x2iSave;
                    }
                }
                else
                {
                    Output("Found factor d=",d," at i=",i); 
                    return d;
                }
            }

            // NOTE: This was not suggested by Pollard's original paper
            k = 0;
            xiSave = xi;
            x2iSave = x2i;
            Qi = 1;
        }
        ++k;
        ++i;
    }
}

vector<BigInt> PollardRhoFactors( const BigInt& n, Int gcdDelay=100 )
{
    vector<BigInt> factors;
    BigInt nRem = n;
    while( true )
    {
        Output("Attempting to factor ",nRem," with a=1");
        PushIndent();
        BigInt factor = PollardRho( nRem, 1, gcdDelay );
        PopIndent();
        if( factor == nRem )
        {
            // Try again with a=-1
            Output("Attempting to factor ",nRem," with a=-1");
            PushIndent();
            factor = PollardRho( nRem, -1, gcdDelay );
            PopIndent();
            factors.push_back( factor );
            if( factor == nRem )
                break; 
            else
                nRem /= factor;
        }
        else
        {
            factors.push_back( factor );
            nRem /= factor;
        }
    }
    return factors;
}
#endif

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

#ifdef EL_HAVE_MPC
    try
    {
        const int minIntBits = Input("--minIntBits","integer bits",256);
        const Int gcdDelay =
          Input("--gcdDelay","GCD delay in Pollard's rho",100);
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        // Try Pollard's rho method with x_0=2 and 
        //   x_{i+1} := x_i^2 - 1 (mod n)
        // (See J.M. Pollard, "A Monte Carlo Method for Factorization")

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) at iterations
        // (100,8200,failed)
        BigInt n = (BigInt(1) << unsigned(77)) - 3;
        Output("n=",n);
        auto factors = PollardRhoFactors( n, gcdDelay );
        Output("Factors of n=",n);
        for( auto factor : factors )
            Output(factor);

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417) at iterations
        // (100,800,5300,failed)
        n = (BigInt(1) << unsigned(79)) - 3;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay );
        Output("Factors of n=",n);
        for( auto factor : factors )
            Output(factor);

        // n = 2^97 - 3
        n = (BigInt(1) << unsigned(97)) - 3;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay );
        Output("Factors of n=",n);
        for( auto factor : factors )
            Output(factor);
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
