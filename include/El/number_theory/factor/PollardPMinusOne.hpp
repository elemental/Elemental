/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_FACTOR_POLLARD_PMINUSONE_HPP
#define EL_NUMBER_THEORY_FACTOR_POLLARD_PMINUSONE_HPP

#ifdef EL_HAVE_MPC
namespace El {

namespace factor {

// Pollard's p-1 can occasionally factor much larger numbers than the rho 
// method but is much less reliable (and the ECM method is a generalization
// which allows it to become more reliable)

namespace pollard_pm1 {

inline BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl )
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
                RuntimeError
                ("Smoothness bound of ",smoothnessBound," exceeded");
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

} // namespace pollard_pm1

inline vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl )
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
        BigInt factor = pollard_pm1::FindFactor( nRem, ctrl );
        PopIndent();
        if( ctrl.time )
            Output("Pollard p-1: ",timer.Stop()," seconds");

        // The factor might be composite, so attempt to factor it
        PushIndent();
        auto subfactors = PollardPMinusOne( factor, ctrl );
        PopIndent();
        for( auto subfactor : subfactors )
            factors.push_back( subfactor );
        nRem /= factor;
    }
    PopIndent();
    sort( factors.begin(), factors.end() );
    return factors;
}

} // namespace factor

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_FACTOR_POLLARDPMINUSONE_HPP
