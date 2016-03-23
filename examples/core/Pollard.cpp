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
                if( d == n )
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
BigInt PollardRhoUnopt( const BigInt& n, Int a=1, Int gcdDelay=100 )
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
                if( d == n )
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

// Aggressively avoid re-allocation
// TODO: Remove the direct mpz_gcd call without sacrificing performance
BigInt PollardRhoOpt( const BigInt& n, Int a=1, Int gcdDelay=100 )
{
    BigInt tmp, d;
    BigInt one(1);

    BigInt xi=2, x2i=2, Qi=1;
    BigInt xiSave=xi, x2iSave=x2i;
    Int k=1, i=1; // it is okay for i to overflow since it is just for printing
    while( true )
    {
        // Advance xi once
        //xi = Mod( xi*xi + a, n );
        tmp = xi;
        tmp *= xi;
        tmp += a;
        xi = tmp;
        xi %= n;

        // Advance x2i twice
        //x2i = Mod( x2i*x2i + a, n );
        tmp = x2i;
        tmp *= x2i;
        tmp += a;
        x2i = tmp;
        x2i %= n;
        //x2i = Mod( x2i*x2i + a, n );
        tmp = x2i;
        tmp *= x2i;
        tmp += a;
        x2i = tmp;
        x2i %= n;

        // Advance Qi
        //Qi = Mod( Qi*(x2i-xi), n );
        tmp = x2i;
        tmp -= xi;
        tmp *= Qi;
        Qi = tmp;
        Qi %= n;

        if( k >= gcdDelay )
        {
            //BigInt d = GCD( Qi, n );
            mpz_gcd( d.Pointer(), Qi.LockedPointer(), n.LockedPointer() );
            if( d > one )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( d == n )
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

BigInt PollardRho( const BigInt& n, Int a=1, Int gcdDelay=100, bool opt=true )
{
    if( opt )
        return PollardRhoOpt( n, a, gcdDelay );
    else
        return PollardRhoUnopt( n, a, gcdDelay );
}

vector<BigInt> PollardRhoFactors
( const BigInt& n,
  Int gcdDelay=100,
  int numReps=20,
  bool opt=true,
  bool time=false )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    Timer timer;
    while( true )
    {
        // Try Miller-Rabin first
        if( time )
            timer.Start();
        Primality primality = PrimalityTest( nRem, numReps );
        if( primality == PRIME )
        {
            if( time )
                Output(nRem," is prime (",timer.Stop()," seconds)");
            else
                Output(nRem," is prime");
            factors.push_back( nRem );       
            break;
        }
        else if( primality == PROBABLY_PRIME )
        {
            if( time )
                Output(nRem," is probably prime (",timer.Stop()," seconds)");
            else
                Output(nRem," is probably prime");
            factors.push_back( nRem );
            break;
        }
        else
        {
            if( time )
                Output(nRem," is composite (",timer.Stop()," seconds)");
            else
                Output(nRem," is composite");
        }

        Output("Attempting to factor ",nRem," with a=1");
        if( time )
            timer.Start();
        PushIndent();
        BigInt factor = PollardRho( nRem, 1, gcdDelay, opt );
        PopIndent();
        if( time )
            Output("Pollard-rho: ",timer.Stop()," seconds");
        if( factor == nRem )
        {
            // Try again with a=-1
            Output("Attempting to factor ",nRem," with a=-1");
            if( time )
                timer.Start();
            PushIndent();
            factor = PollardRho( nRem, -1, gcdDelay, opt );
            PopIndent();
            if( time )
                Output("Pollard-rho: ",timer.Stop()," seconds");
            // TODO: Test for (unlikely) possibility of composite factor?
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

BigInt PollardPMinusOneUnopt( const BigInt& n, BigInt smoothness=100000 )
{
    double nLog = Log( n );
    BigInt smoothnessBound = Pow( BigInt(10), BigInt(9) );
    while( true )
    {
        // Uniformly select a in (Z/(n))*
        // (alternatively, we could set a=2)
        BigInt a = SampleUniform( BigInt(0), n );
        while( GCD( a, n ) != BigInt(1) )
        {
            a = SampleUniform( BigInt(0), n ); 
        }
        Output("a=",a);

        // TODO: Pre-sieve the list of primes rather than using NextPrime

        // TODO: Enable this until a gcd of n is found
        /*
        BigInt q(2), b(a);
        while( q <= smoothness )
        {
            unsigned qExponent = unsigned(nLog/Log(q));
            for( Int i=0; i<qExponent; ++i )
                b = PowMod( b, q, n );
            q = NextPrime( q );
        }
        BigInt g = GCD( b-1, n );
        */

        // Divide by 2 last
        BigInt q(3), b(a);
        while( q <= smoothness )
        {
            unsigned qExponent = unsigned(nLog/Log(q));
            for( Int i=0; i<qExponent; ++i )
                b = PowMod( b, q, n );
            q = NextPrime( q );
        }
        unsigned twoExponent = unsigned(nLog/Log(2));
        BigInt g = GCD( b-1, n );
        for( Int i=0; i<twoExponent; ++i )
        {
            if( g > BigInt(1) && g < n )
            {
                Output("Found factor ",g," at i=",i);
                break;
            }
            b = PowMod( b, BigInt(2), n );
            g = GCD( b-1, n );
        }

        // TODO: Add support for a second stage

        if( g == BigInt(1) )
        {
            if( smoothness >= smoothnessBound )
            {
                Output("Smoothness bound of ",smoothnessBound," exceeded");
                return n;
            }
            smoothness *= 2;
            Output("Increased smoothness to ",smoothness);
        }
        else if( g == n )
        {
            if( smoothness <= BigInt(2) )
            {
                Output("Smoothness lower bound of 2 hit");
                return n;
            }
            smoothness /= 2;
            Output("Decreased smoothness to ",smoothness);
        }
        else
        {
            Output("Found factor of ",g);
            return g;
        }
    }
}

BigInt PollardPMinusOneOpt( const BigInt& n, BigInt smoothness=100000 )
{
    double nLog = Log( n );
    BigInt smoothnessBound = Pow( BigInt(10), BigInt(9) );
    while( true )
    {
        // Uniformly select a in (Z/(n))*
        // (alternatively, we could set a=2)
        BigInt a = SampleUniform( BigInt(0), n );
        while( GCD( a, n ) != BigInt(1) )
        {
            a = SampleUniform( BigInt(0), n ); 
        }
        Output("a=",a);

        /*
        BigInt q(2), b(a);
        while( q <= smoothness )
        {
            unsigned qExponent = unsigned(nLog/Log(q));
            for( Int i=0; i<qExponent; ++i )
                b = PowMod( b, q, n );
            q = NextPrime( q );
        }
        BigInt g = GCD( b-1, n );
        */

        // Divide by 2 last
        BigInt q(3), b(a);
        while( q <= smoothness )
        {
            unsigned qExponent = unsigned(nLog/Log(q));
            for( Int i=0; i<qExponent; ++i )
            {
                //b = PowMod( b, q, n );
                mpz_powm
                ( b.Pointer(),
                  b.LockedPointer(),
                  q.LockedPointer(),
                  n.LockedPointer() );
            }
            //q = NextPrime( q );
            mpz_nextprime( q.Pointer(), q.LockedPointer() );
        }
        unsigned twoExponent = unsigned(nLog/Log(2));
        BigInt g = GCD( b-1, n );
        BigInt two(2), tmp;
        for( Int i=0; i<twoExponent; ++i )
        {
            if( g > BigInt(1) && g < n )
            {
                Output("Found factor ",g," at i=",i);
                break;
            }
            //b = PowMod( b, two, n );
            mpz_powm
            ( b.Pointer(),
              b.LockedPointer(),
              two.LockedPointer(),
              n.LockedPointer() );
            //g = GCD( b-1, n );
            tmp = b;
            tmp -= 1;
            mpz_gcd( g.Pointer(), tmp.LockedPointer(), n.LockedPointer() );
        }

        if( g == BigInt(1) )
        {
            if( smoothness >= smoothnessBound )
            {
                Output("Smoothness bound of ",smoothnessBound," exceeded");
                return n;
            }
            smoothness *= 2;
            Output("Increased smoothness to ",smoothness);
        }
        else if( g == n )
        {
            if( smoothness <= BigInt(2) )
            {
                Output("Smoothness lower bound of 2 hit");
                return n;
            }
            smoothness /= 2;
            Output("Decreased smoothness to ",smoothness);
        }
        else
        {
            Output("Found factor of ",g);
            return g;
        }
    }
}

BigInt PollardPMinusOne
( const BigInt& n, BigInt smoothness=100000, bool opt=true )
{
    if( opt )
        return PollardPMinusOneOpt( n, smoothness );
    else
        return PollardPMinusOneUnopt( n, smoothness );
}

vector<BigInt> PollardPMinusOneFactors
( const BigInt& n,
  BigInt smoothness=100000,
  int numReps=20,
  bool opt=true,
  bool time=false )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    Timer timer;
    while( true )
    {
        // Try Miller-Rabin first
        if( time )
            timer.Start();
        Primality primality = PrimalityTest( nRem, numReps );
        if( primality == PRIME )
        {
            if( time )
                Output(nRem," is prime (",timer.Stop()," seconds)");
            else
                Output(nRem," is prime");
            factors.push_back( nRem );       
            break;
        }
        else if( primality == PROBABLY_PRIME )
        {
            if( time )
                Output(nRem," is probably prime (",timer.Stop()," seconds)");
            else
                Output(nRem," is probably prime");
            factors.push_back( nRem );
            break;
        }
        else
        {
            if( time )
                Output(nRem," is composite (",timer.Stop()," seconds)");
            else
                Output(nRem," is composite");
        }

        Output("Attempting to factor ",nRem);
        if( time )
            timer.Start();
        PushIndent();
        BigInt factor = PollardPMinusOne( nRem, smoothness, opt );
        PopIndent();
        if( time )
            Output("Pollard p-1: ",timer.Stop()," seconds");

        if( factor == nRem )
        {
            factors.push_back( factor );
            break; 
        }
        else
        {
            // The factor might be composite, so attempt to factor it
            PushIndent();
            auto subfactors =
              PollardPMinusOneFactors( factor, smoothness, numReps, opt, time );
            PopIndent();
            for( auto subfactor : subfactors )
                factors.push_back( subfactor );
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
        const int minIntBits = Input("--minIntBits","integer bits",384);
        const Int gcdDelay =
          Input("--gcdDelay","GCD delay in Pollard's rho",100);
        const int numReps = Input("--numReps","num Miller-Rabin reps,",20);
        const bool opt = Input("--opt","optimized allocations?",true);
        const bool time = Input("--time","time Pollard rho steps?",true);
        const BigInt smoothness =
          Input("--smoothness","smoothness bound for p-1",BigInt(100000));
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        // Try Pollard's rho method with x_0=2 and 
        //   x_{i+1} := x_i^2 - 1 (mod n)
        // (See J.M. Pollard, "A Monte Carlo Method for Factorization")

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) at iterations
        // (100,8200,failed)
        BigInt n = Pow(BigInt(2),unsigned(77)) - 3;
        Output("n=",n);
        auto factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");
        factors = PollardPMinusOneFactors( n, smoothness, numReps, opt, time );
        Output("");

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417) at iterations
        // (100,800,5300,failed)
        n = Pow(BigInt(2),unsigned(79)) - 3;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");
        factors = PollardPMinusOneFactors( n, smoothness, numReps, opt, time );
        Output("");

        // n = 2^97 - 3
        n = Pow(BigInt(2),unsigned(97)) - 3;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");
        factors = PollardPMinusOneFactors( n, smoothness, numReps, opt, time );
        Output("");

        // n = 3^100 + 2 (example from http://julia-programming-language.2336112.n4.nabble.com/Factorization-of-big-integers-is-taking-too-long-td15925.html)
        n = Pow(BigInt(3),unsigned(100)) + 2;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");
        factors = PollardPMinusOneFactors( n, smoothness, numReps, opt, time );
        Output("");

        // TODO: Re-enable this when the methods are more refined
        // n = 3^100 + 4
        /* 
        n = Pow(BigInt(3),unsigned(100)) + 4;
        Output("n=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");
        factors = PollardPMinusOneFactors( n, smoothness, numReps, opt, time );
        */
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
