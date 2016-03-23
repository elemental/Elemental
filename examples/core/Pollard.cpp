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
            BigInt gcd = GCD( Qi, n );
            if( gcd > BigInt(1) )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( gcd == n )
                {
                    if( gcdDelay == 1 )
                    {
                        Output("(x_n) converged before (x_n mod p) at i=",i);
                        return gcd;
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
                    Output("Found factor ",gcd," at i=",i); 
                    return gcd;
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

// TODO: Remove the direct mpz_gcd call without sacrificing performance
BigInt PollardRhoOpt( const BigInt& n, Int a=1, Int gcdDelay=100 )
{
    BigInt tmp, gcd;
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
            //BigInt gcd = GCD( Qi, n );
            mpz_gcd( gcd.Pointer(), Qi.LockedPointer(), n.LockedPointer() );
            if( gcd > one )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( gcd == n )
                {
                    if( gcdDelay == 1 )
                    {
                        Output("(x_n) converged before (x_n mod p) at i=",i);
                        return gcd;
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
                    Output("Found factor ",gcd," at i=",i); 
                    return gcd;
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

// TODO: Unoptimized version

BigInt PollardPMinusOneOpt( const BigInt& n, BigInt smoothness=100000 )
{
    const double twoLog = Log( 2 );
    const double nLog = Log( n );
    const BigInt zero(0), one(1);
    const BigInt smoothnessBound =
        Max( Pow(BigInt(10),BigInt(9)), smoothness*8 );

    bool separateOdd=false;

    BigInt smallPrime, gcd, tmp;
    while( true )
    {
        // Uniformly select a in (Z/(n))*
        // (alternatively, we could set a=2)
        BigInt a = SampleUniform( zero, n );
        while( GCD( a, n ) != one )
        {
            a = SampleUniform( zero, n ); 
        }

        // TODO: Generate the list of primes with a sieve instead of nextprime
        smallPrime = ( separateOdd ? 3 : 2 );
        while( smallPrime <= smoothness )
        {
            unsigned smallPrimeExponent = unsigned(nLog/Log(smallPrime));
            for( Int i=0; i<smallPrimeExponent; ++i )
            {
                //a = PowMod( a, smallPrime, n );
                mpz_powm
                ( a.Pointer(),
                  a.LockedPointer(),
                  smallPrime.LockedPointer(),
                  n.LockedPointer() );
            }
            //smallPrime = NextPrime( smallPrime );
            mpz_nextprime( smallPrime.Pointer(), smallPrime.LockedPointer() );
        }

        // gcd := GCD( a-1, n )
        tmp = a; 
        tmp -= 1;
        mpz_gcd( gcd.Pointer(), tmp.LockedPointer(), n.LockedPointer() );

        if( separateOdd )
        { 
            unsigned twoExponent = unsigned(nLog/twoLog);
            for( Int i=0; i<twoExponent; ++i )
            {
                if( gcd > one && gcd < n )
                    break;
                //a = PowMod( a, 2, n );
                a *= a;
                a %= n;
                //gcd = GCD( a-1, n );
                tmp = a;
                tmp -= 1;
                mpz_gcd( gcd.Pointer(), tmp.LockedPointer(), n.LockedPointer() );
            }
        }

        if( gcd == one )
        {
            if( smoothness >= smoothnessBound )
            {
                Output("Smoothness bound of ",smoothnessBound," exceeded");
                return n;
            }
            smoothness *= 2;
            Output("Increased smoothness to ",smoothness);
        }
        else if( gcd == n )
        {
            separateOdd = true;
            Output("Separately checking powers of two");
        }
        else
        {
            Output("Found factor of ",gcd);
            return gcd;
        }
    }
}

// Pollard's p-1 can occasionally factor much larger numbers than the rho 
// method but is much less reliable (and the ECM method is a generalization
// which allows it to become more reliable)
BigInt PollardPMinusOne
( const BigInt& n, BigInt smoothness=100000, bool opt=true )
{
    return PollardPMinusOneOpt( n, smoothness );
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
        const bool largeRho =
          Input("--largeRho","reproduce Pollard and Brent's result?",false);
        const bool heroicPMinusOne =
          Input("--heroicPMinusOne","reproduce Zimmerman's result?",false);
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) 
        BigInt n = Pow(BigInt(2),unsigned(77)) - 3;
        Output("n=2^77-3=",n);
        auto factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417)
        n = Pow(BigInt(2),unsigned(79)) - 3;
        Output("n=2^79-3=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");

        // n = 2^97 - 3
        n = Pow(BigInt(2),unsigned(97)) - 3;
        Output("n=2^97-3=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
        Output("");

        // n = 3^100 + 2 
        n = Pow(BigInt(3),unsigned(100)) + 2;
        Output("n=3^100+2=",n);
        factors = PollardRhoFactors( n, gcdDelay, numReps, opt, time );
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
            factors = PollardRhoFactors( z, gcdDelay, numReps, opt, time );
            Output("");
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
            BigInt smoothness = Pow(BigInt(10),BigInt(13));
            Output("z=2^2098+1=",z);
            factors =
              PollardPMinusOneFactors( z, smoothness, numReps, opt, time );
        }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
