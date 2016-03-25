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
struct PollardRhoCtrl
{
    Int a0=1;
    Int a1=-1;
    unsigned long numSteps=1u;
    BigInt x0=BigInt(2);
    Int gcdDelay=100;
    int numReps=30;
    bool progress=false;
    bool time=false;
};

// TODO: Add the ability to set a maximum number of iterations
BigInt PollardRhoInner
( const BigInt& n,
  Int a=1,
  const PollardRhoCtrl& ctrl=PollardRhoCtrl() )
{
    if( a == 0 || a == -2 )
        Output("WARNING: Problematic choice of Pollard rho shift");
    BigInt tmp, gcd;
    BigInt one(1);

    auto xAdvance =
      [&]( BigInt& x )
      {
        if( ctrl.numSteps == 1 )
        {
            // TODO: Determine if there is a penalty to x *= x
            /*
            tmp = x;
            tmp *= x;
            tmp += a;
            x = tmp;
            x %= n;
            */
            x *= x;
            x += a;
            x %= n;
        }
        else
        {
            PowMod( x, 2*ctrl.numSteps, n, x );
            x += a;
            x %= n;
        }
      };

    auto QAdvance =
      [&]( const BigInt& x, const BigInt& x2, BigInt& Q )
      {
        tmp = x2;
        tmp -= x;
        Q *= tmp;
        Q %= n;
      };

    Int gcdDelay = ctrl.gcdDelay;
    BigInt xi=ctrl.x0, x2i=ctrl.x0, Qi=1;
    BigInt xiSave=xi, x2iSave=x2i;
    Int k=1, i=1; // it is okay for i to overflow since it is just for printing
    while( true )
    {
        // Advance xi once
        xAdvance( xi );

        // Advance x2i twice
        xAdvance( x2i );
        xAdvance( x2i );

        // Advance Qi
        QAdvance( xi, x2i, Qi );

        if( k >= gcdDelay )
        {
            GCD( Qi, n, gcd );
            if( gcd > one )
            {
                // NOTE: This was not suggested by Pollard's original paper
                if( gcd == n )
                {
                    if( gcdDelay == 1 )
                    {
                        if( ctrl.progress )
                            Output("(x) converged before (x mod p) at i=",i);
                        return gcd;
                    }
                    else
                    {
                        if( ctrl.progress )
                            Output("Backtracking at i=",i);
                        i = Max( i-(gcdDelay+1), 0 );
                        gcdDelay = 1;
                        xi = xiSave;
                        x2i = x2iSave;
                    }
                }
                else
                {
                    if( ctrl.progress )
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

vector<BigInt> PollardRho
( const BigInt& n,
  const PollardRhoCtrl& ctrl=PollardRhoCtrl() )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    Timer timer;
    PushIndent();
    while( true )
    {
        // Try Miller-Rabin first
        if( ctrl.time )
            timer.Start();
        Primality primality = PrimalityTest( nRem, ctrl.numReps );
        if( primality == PRIME )
        {
            if( ctrl.time )
                Output(nRem," is prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is prime");
            factors.push_back( nRem );       
            break;
        }
        else if( primality == PROBABLY_PRIME )
        {
            if( ctrl.time )
                Output(nRem," is probably prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is probably prime");
            factors.push_back( nRem );
            break;
        }
        else
        {
            if( ctrl.time )
                Output(nRem," is composite (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is composite");
        }

        if( ctrl.progress )
            Output("Attempting to factor ",nRem," with a=",ctrl.a0);
        if( ctrl.time )
            timer.Start();
        PushIndent();
        BigInt factor = PollardRhoInner( nRem, ctrl.a0, ctrl );
        PopIndent();
        if( ctrl.time )
            Output("Pollard-rho (a=",ctrl.a0,"): ",timer.Stop()," seconds");
        if( factor == nRem )
        {
            // Try again with a=ctrl.a1
            if( ctrl.progress )
                Output("Attempting to factor ",nRem," with a=",ctrl.a1);
            if( ctrl.time )
                timer.Start();
            PushIndent();
            factor = PollardRhoInner( nRem, ctrl.a1, ctrl );
            PopIndent();
            if( ctrl.time )
                Output("Pollard-rho (a=",ctrl.a1,"): ",timer.Stop()," seconds");
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
    PopIndent();
    sort( factors.begin(), factors.end() );
    return factors;
}

struct PollardPMinusOneCtrl
{
    BigInt smoothness=BigInt(100000);
    int numReps=30;
    bool progress=false;
    bool time=false;
};

// Pollard's p-1 can occasionally factor much larger numbers than the rho 
// method but is much less reliable (and the ECM method is a generalization
// which allows it to become more reliable)
BigInt PollardPMinusOneInner
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() )
{
    const double twoLog = Log( 2 );
    const double nLog = Log( n );
    const BigInt zero(0), one(1);

    BigInt smoothness = ctrl.smoothness;
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
                // a = a^smallPrime (mod n)
                PowMod( a, smallPrime, n, a );
            }
            // Move smallPrime to the next higher prime
            NextPrime( smallPrime, smallPrime );
        }

        // gcd := GCD( a-1, n )
        tmp = a; 
        tmp -= 1;
        GCD( tmp, n, gcd );

        if( separateOdd )
        { 
            unsigned twoExponent = unsigned(nLog/twoLog);
            for( Int i=0; i<twoExponent; ++i )
            {
                if( gcd > one && gcd < n )
                    break;
                // a = a*a (mod n)
                a *= a;
                a %= n;
                //gcd = GCD( a-1, n );
                tmp = a;
                tmp -= 1;
                GCD( tmp, n, gcd );
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
            if( ctrl.progress )
                Output("Increased smoothness to ",smoothness);
        }
        else if( gcd == n )
        {
            separateOdd = true;
            if( ctrl.progress )
                Output("Separately checking powers of two");
        }
        else
        {
            if( ctrl.progress )
                Output("Found factor of ",gcd);
            return gcd;
        }
    }
}

vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    Timer timer;
    PushIndent();
    while( true )
    {
        // Try Miller-Rabin first
        if( ctrl.time )
            timer.Start();
        Primality primality = PrimalityTest( nRem, ctrl.numReps );
        if( primality == PRIME )
        {
            if( ctrl.time )
                Output(nRem," is prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is prime");
            factors.push_back( nRem );       
            break;
        }
        else if( primality == PROBABLY_PRIME )
        {
            if( ctrl.time )
                Output(nRem," is probably prime (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is probably prime");
            factors.push_back( nRem );
            break;
        }
        else
        {
            if( ctrl.time )
                Output(nRem," is composite (",timer.Stop()," seconds)");
            else if( ctrl.progress )
                Output(nRem," is composite");
        }

        if( ctrl.progress )
            Output("Attempting to factor ",nRem);
        if( ctrl.time )
            timer.Start();
        PushIndent();
        BigInt factor = PollardPMinusOneInner( nRem, ctrl );
        PopIndent();
        if( ctrl.time )
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
            auto subfactors = PollardPMinusOne( factor, ctrl );
            PopIndent();
            for( auto subfactor : subfactors )
                factors.push_back( subfactor );
            nRem /= factor;
        }
    }
    PopIndent();
    sort( factors.begin(), factors.end() );
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
        const unsigned long numSteps =
          Input("--numSteps","x^2k + a in Pollard rho",1u);
        const BigInt x0 = Input("--x0","x0 in Pollard rho",BigInt(2));
        const Int gcdDelay =
          Input("--gcdDelay","GCD delay in Pollard's rho",100);
        const int numReps = Input("--numReps","num Miller-Rabin reps,",30);
        const bool progress = Input("--progress","factor progress?",true);
        const bool time = Input("--time","time Pollard rho steps?",true);
        const bool largeRho =
          Input("--largeRho","reproduce Pollard and Brent's result?",false);
        const bool heroicPMinusOne =
          Input("--heroicPMinusOne","reproduce Zimmerman's result?",false);
        ProcessInput();
        PrintInputReport(); 

        mpc::SetMinIntBits( minIntBits );

        PollardRhoCtrl rhoCtrl;
        rhoCtrl.numSteps = numSteps;
        rhoCtrl.x0 = x0;
        rhoCtrl.gcdDelay = gcdDelay;
        rhoCtrl.numReps = numReps;
        rhoCtrl.progress = progress;
        rhoCtrl.time = time;

        PollardPMinusOneCtrl pm1Ctrl;
        pm1Ctrl.numReps = numReps;
        pm1Ctrl.progress = progress;
        pm1Ctrl.time = time;

        // n = 2^77 - 3
        // We should find (1291,99432527,1177212722617) 
        BigInt n = Pow(BigInt(2),unsigned(77)) - 3;
        Output("n=2^77-3=",n);
        auto factors = PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 2^79 - 3
        // We should find (5,3414023,146481287,241741417)
        n = Pow(BigInt(2),unsigned(79)) - 3;
        Output("n=2^79-3=",n);
        factors = PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 2^97 - 3
        n = Pow(BigInt(2),unsigned(97)) - 3;
        Output("n=2^97-3=",n);
        factors = PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
        Output("");

        // n = 3^100 + 2 
        n = Pow(BigInt(3),unsigned(100)) + 2;
        Output("n=3^100+2=",n);
        factors = PollardRho( n, rhoCtrl );
        Output("factors:");
        for( auto factor : factors )
            Output("  ",factor); 
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
            factors = PollardRho( z, rhoCtrl );
            Output("factors:");
            for( auto factor : factors )
                Output("  ",factor); 
            Output("");

            if( numSteps != 512u )
            {
                // Try again with numSteps=512 and x0=3
                Output("Factoring z=2^(2^8)+1 again with numSteps=512");
                auto rhoCtrlMod = rhoCtrl;
                rhoCtrlMod.numSteps = 512u;
                rhoCtrlMod.x0 = BigInt(3);
                factors = PollardRho( z, rhoCtrlMod );
                Output("factors:");
                for( auto factor : factors )
                    Output("  ",factor); 
                Output("");
            }
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
            pm1Ctrl.smoothness = Pow(BigInt(10),BigInt(13));
            Output("z=2^2098+1=",z);
            factors = PollardPMinusOne( z, pm1Ctrl );
            Output("factors:");
            for( auto factor : factors )
                Output("  ",factor); 
        }
    }
    catch( std::exception& e ) { ReportException(e); }
#endif

    return 0;
}
