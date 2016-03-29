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
        DynamicSieve<unsigned long long>& sieve,
  const PollardPMinusOneCtrl& ctrl )
{
    const double twoLog = Log( 2. );
    const double nLog = double( Log( BigFloat(n) ) );
    const BigInt zero(0), one(1);

    unsigned long long smooth1 = ctrl.smooth1;
    unsigned long long smooth2 = ctrl.smooth2;
    const auto smooth1Bound = Max( Pow(10ULL,9ULL), 8ULL*smooth1 );
    const auto smooth2Bound = Max( Pow(10ULL,10ULL), 8ULL*smooth2 );

    // Ensure that we have sieved at least up until smooth1
    while( sieve.oddPrimes.back() < smooth1 )
        sieve.NextPrime();

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

        if( !separateOdd )
        {
            // Handle 2 separately
            unsigned smallPrimeExponent = unsigned(nLog/twoLog);
            for( Int i=0; i<smallPrimeExponent; ++i )
            {
                // a = a^2 (mod n)
                a *= a;
                a %= n;
            }
        }
        auto smooth1Iter = 
          std::upper_bound
          ( sieve.oddPrimes.begin(), sieve.oddPrimes.end(), smooth1 );
        for( auto iter=sieve.oddPrimes.begin(); iter<smooth1Iter; ++iter )
        {
            auto smallPrime = *iter;
            double smallPrimeLog = double(Log(double(smallPrime)));
            unsigned smallPrimeExponent = unsigned(nLog/smallPrimeLog);
            for( Int i=0; i<smallPrimeExponent; ++i )
            {
                // a = a^smallPrime (mod n)
                PowMod( a, smallPrime, n, a );
            }
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

        // TODO: Introduce second stage

        if( gcd == one )
        {
            if( smooth1 >= smooth1Bound )
            {
                RuntimeError
                ("Stage-1 smoothness bound of ",smooth1Bound," exceeded");
            }
            smooth1 *= 2;
            if( ctrl.progress )
                Output("Increased stage-1 smoothness to ",smooth1);
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

inline BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl )
{
    DynamicSieve<unsigned long long> sieve;
    return FindFactor( n, sieve, ctrl );
}

} // namespace pollard_pm1

inline vector<BigInt> PollardPMinusOne
( const BigInt& n,
        DynamicSieve<unsigned long long>& sieve,
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
        BigInt factor = pollard_pm1::FindFactor( nRem, sieve, ctrl );
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

inline vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl )
{
    DynamicSieve<unsigned long long> sieve;
    return PollardPMinusOne( n, sieve, ctrl );
}

} // namespace factor

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_FACTOR_POLLARDPMINUSONE_HPP
