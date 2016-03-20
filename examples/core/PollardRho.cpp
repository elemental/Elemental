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
// TODO: Generalize beyond x_{i+1} = x_i^2 - 1 (mod n)
BigInt PollardRho( const BigInt& n, Int m, Int S )
{
    BigInt xi=2, x2i=2, Qi=1;
    for( Int i=0; i<S; ++i )
    {
        // Advance xi once
        xi = Mod( xi*xi - 1, n );

        // Advance x2i twice
        x2i = Mod( x2i*x2i - 1, n );
        x2i = Mod( x2i*x2i - 1, n );

        // Advance Qi
        Qi = Mod( Qi*(x2i-xi), n );

        if( i % m == 0 )
        {
            BigInt d = GCD( Qi, n );
            if( BigInt(1) < d && d < n )
            {
                Output("Found factor d=",d," at i=",i); 
                return d;
            }
        }
        if( i == S-1 )
        {
            Output("Failed to find a factor in time");
            Output("Final xi=",xi,", x2i=",x2i,", and Qi=",Qi);
        }
    }
    return n;
}

vector<BigInt> PollardRhoFactors( const BigInt& n, Int m, Int S )
{
    vector<BigInt> factors;
    BigInt a = n;
    while( true )
    {
        BigInt factor = PollardRho( a, m, S );
        factors.push_back( factor );
        if( factor == a )
            break; 
        else
            a /= factor;
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
        const int minIntBits = Input("--minIntBits","integer bits",128);
        const Int m = Input("--m","m from Pollard's rho",100);
        const Int S = Input("--S","S from Pollard's rho",10000);
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        // Try Pollard's rho method with x_0=2 and 
        //   x_{i+1} := x_i^2 - 1 (mod n)
        // (See J.M. Pollard, "A Monte Carlo Method for Factorization")

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) at iterations
        // (100,8200,failed)
        BigInt n( 1 );
        n <<= unsigned(77);
        n -= 3;
        Output("n=",n);
        auto factors = PollardRhoFactors( n, m, S );
        Output("Factors of n=",n);
        for( auto factor : factors )
            Output(factor);

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417) at iterations
        // (100,800,5300,failed)
        n = 1;
        n <<= unsigned(79);
        n -= 3;
        Output("n=",n);
        factors = PollardRhoFactors( n, m, S );
        Output("Factors of n=",n);
        for( auto factor : factors )
            Output(factor);
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
