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

template<typename TSieve,typename TSieveSmall>
inline BigInt FindFactor
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    const double twoLog = Log( 2. );
    const double nLog = double( Log( BigFloat(n) ) );
    const BigInt zero(0), one(1);

    TSieve smooth1 = ctrl.smooth1;
    TSieve smooth2 = ctrl.smooth2;
    TSieve smooth1Bound = Max( Pow(10ULL,9ULL), 8ULL*smooth1 );
    TSieve smooth2Bound = Max( Pow(10ULL,10ULL), 8ULL*smooth2 );

    // Ensure that we do not have a GCD of n appear too many times despite
    // separately checking the powers of two
    bool separateOdd=false;
    Int maxGCDFailures=10;
    Int numGCDFailures=0;

    BigInt smallPrime, gcd, tmp, diffPower;
    while( true )
    {
        // Ensure that we have sieved at least up until smooth1
        bool neededStage1Sieving = ( sieve.oddPrimes.back() < smooth1 );
        if( ctrl.progress && neededStage1Sieving )
            Output
            ("Updating sieve from ",sieve.oddPrimes.back()," to ",smooth1);
        sieve.Generate( smooth1 );
        if( ctrl.progress && neededStage1Sieving )
            Output("Done sieving for stage 1");

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
        auto smooth1End = 
          std::upper_bound
          ( sieve.oddPrimes.begin(),
            sieve.oddPrimes.end(),
            smooth1 );
        for( auto iter=sieve.oddPrimes.begin(); iter<smooth1End; ++iter )
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
        if( gcd > one && gcd < n )
        {
            if( ctrl.progress )
                Output("Found stage-1 factor of ",gcd);
            return gcd;
        }

        if( separateOdd )
        { 
            unsigned twoExponent = unsigned(nLog/twoLog);
            for( Int i=0; i<twoExponent; ++i )
            {
                // a = a*a (mod n)
                a *= a;
                a %= n;

                //gcd = GCD( a-1, n );
                tmp = a;
                tmp -= 1;
                GCD( tmp, n, gcd );
                if( gcd > one && gcd < n )
                {
                    if( ctrl.progress )
                        Output("Found separate stage-1 factor of ",gcd);
                    return gcd;
                }
            }
        }

        if( gcd == n )
        {
            if( separateOdd )
            {
                ++numGCDFailures;
                if( numGCDFailures >= maxGCDFailures )
                    RuntimeError
                    ("Too many GCD failures despite separately checking "
                     "powers of two");
            }
            else
            {
                ++numGCDFailures;
                separateOdd = true;
                if( ctrl.progress )
                    Output("GCD was n; will separately check powers of two");
            }
            continue;
        }
        else // gcd == one
        {
            if( ctrl.progress )
                Output("Proceeding to stage-2");
        }
        
        // Run stage-2
        // ----------- 

        // Store all powers of a^{d_i}, where d_i is the difference between
        // primes p_{i+1} and p_i, where smooth1 < p_{i+1} <= smooth2
        bool neededStage2Sieving = ( sieve.oddPrimes.back() < smooth2 );
        if( ctrl.progress && neededStage2Sieving )
            Output
            ("Updating sieve from ",sieve.oddPrimes.back()," to ",smooth2);
        sieve.Generate( smooth2 );
        if( ctrl.progress && neededStage2Sieving )
            Output("Done sieving for stage 2");

        // NOTE: stage1End has potentially been invalidated due to reallocation
        auto stage2Beg =
          std::lower_bound
          ( sieve.oddPrimes.begin(),
            sieve.oddPrimes.end(),
            smooth1 );
        auto stage2End =
          std::upper_bound
          ( stage2Beg,
            sieve.oddPrimes.end(),
            smooth2 );
        std::map<TSieve,BigInt> diffPowers;
        for( auto iter=stage2Beg; iter<stage2End; ++iter )
        {
            TSieve diff = *iter - *(iter-1);
            auto search = diffPowers.find( diff );

            const size_t whichPrime = (iter-sieve.oddPrimes.begin())+1;
            Output("  ",whichPrime,": ",*iter," - ",*(iter-1)," = ",diff);
            if( search == diffPowers.end() )
            {
                PowMod( a, diff, n, diffPower );
                diffPowers.insert( std::make_pair(diff,diffPower) );
                Output("    Stored ",a,"^",diff,"=",diffPower);
            }
            else
                Output("    diff=",diff," was redundant");
        }

        // Test each stage two candidate
        for( auto iter=stage2Beg; iter<stage2End; ++iter )
        {
            TSieve diff = *iter - *(iter-1);
            a *= diffPowers[diff];
            a %= n;
            
            //gcd = GCD( a-1, n );
            tmp = a;
            tmp -= 1;
            GCD( tmp, n, gcd );
            if( gcd > one && gcd < n )
            {
                if( ctrl.progress )
                    Output("Found stage-2 factor of ",gcd);
                return gcd;
            }
        }

        if( gcd == n )
        {
            if( separateOdd )
            {
                if( ctrl.progress )
                    Output("GCD failure ",numGCDFailures);
                ++numGCDFailures;
                if( numGCDFailures >= maxGCDFailures )
                    RuntimeError
                    ("Too many GCD failures despite separately checking "
                     "powers of two");
            }
            else
            {
                ++numGCDFailures;
                separateOdd = true;
                if( ctrl.progress )
                    Output("GCD was n; will separately check powers of two");
            }
            continue;
        }
        else // gcd == one
        {
            if( smooth1 >= smooth1Bound )
            {
                RuntimeError
                ("Stage-1 smoothness bound of ",smooth1Bound," exceeded");
            }
            if( smooth2 >= smooth2Bound )
            {
                RuntimeError
                ("Stage-2 smoothness bound of ",smooth2Bound," exceeded");
            }
            smooth1 *= 2;
            smooth2 *= 2;
            if( ctrl.progress )
                Output
                ("Increased stage-1 smoothness to ",smooth1," and stage-2 "
                 "smoothness to ",smooth2);
        }
    }
}

template<typename TSieve,typename TSieveSmall>
inline BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    DynamicSieve<TSieve,TSieveSmall> sieve;
    return FindFactor( n, sieve, ctrl );
}

} // namespace pollard_pm1

template<typename TSieve,typename TSieveSmall>
inline vector<BigInt> PollardPMinusOne
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
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

template<typename TSieve,typename TSieveSmall>
inline vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    DynamicSieve<TSieve,TSieveSmall> sieve;
    return PollardPMinusOne( n, sieve, ctrl );
}

} // namespace factor

} // namespace El

#endif // ifdef EL_HAVE_MPC

#endif // ifndef EL_NUMBER_THEORY_FACTOR_POLLARDPMINUSONE_HPP
