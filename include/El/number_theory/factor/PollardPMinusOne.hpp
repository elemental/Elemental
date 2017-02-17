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

inline void RepeatedSquareMod
( BigInt& a,
  const BigInt& n,
  const double& nLog )
{
    double twoLog = Log(2.);
    unsigned exponent = unsigned(nLog/twoLog);
    for( unsigned i=0; i<exponent; ++i )
    {
        a *= a;
        a %= n;
    }
}

template<typename TSieve>
void RepeatedPowMod
( BigInt& a,
  TSieve p,
  const BigInt& n,
  const double& nLog )
{
    double pLog = double(Log(double(p)));
    unsigned exponent = unsigned(nLog/pLog);
    for( unsigned i=0; i<exponent; ++i )
        PowMod( a, p, n, a );
}

template<typename Iterator>
void RepeatedPowModRange
( BigInt& a,
  Iterator pBeg,
  Iterator pEnd,
  const BigInt& n,
  const double& nLog,
  bool checkpoint=false,
  Int checkpointFreq=1000000 )
{
    Int checkpointCounter = 0;
    for( auto pPtr=pBeg; pPtr<pEnd; ++pPtr )
    {
        auto p = *pPtr;
        RepeatedPowMod( a, p, n, nLog );

        ++checkpointCounter;
        if( checkpoint && checkpointCounter >= checkpointFreq )
        {
            Output("After p=",p,", exponential was a=",a);
            checkpointCounter = 0;
        }
    }
}

// NOTE: Returns the GCD of stage 1 and overwrites a with a power of a
template<typename TSieve,typename TSieveSmall>
BigInt StageOne
( const BigInt& n,
        BigInt& a,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
        bool separateOdd,
        TSieve primeBound,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    const BigInt& one = BigIntOne();

    // Ensure that we have sieved at least up until primeBound
    bool neededStage1Sieving = ( sieve.oddPrimes.back() < primeBound );
    if( ctrl.progress )
    {
        if( neededStage1Sieving )
            Output
            ("Updating sieve from ",sieve.oddPrimes.back()," to ",primeBound);
        else
            Output("Stage 1 sieve was already sufficient");
    }
    sieve.Generate( primeBound );
    if( ctrl.progress && neededStage1Sieving )
        Output("Done sieving for stage 1");

    double nLog = double(Log(BigFloat(n)));
    if( !separateOdd && !ctrl.jumpstart1 )
    {
        RepeatedSquareMod( a, n, nLog );
    }
    auto oddPrimeBeg = sieve.oddPrimes.begin();
    if( ctrl.jumpstart1 )
    {
        // Start at the first odd prime >= ctrl.start1
        oddPrimeBeg =
          std::lower_bound
          ( oddPrimeBeg, sieve.oddPrimes.end(), ctrl.start1 );
    }
    // We stop just before the last prime *greater* than primeBound
    auto oddPrimeEnd =
      std::upper_bound( oddPrimeBeg, sieve.oddPrimes.end(), primeBound );
    RepeatedPowModRange
    ( a, oddPrimeBeg, oddPrimeEnd, n, nLog,
      ctrl.checkpoint, ctrl.checkpointFreq );
    if( ctrl.progress )
        Output("Done with stage-1 exponentiation");

    // gcd := GCD( a-1, n )
    BigInt tmp = a;
    tmp -= 1;
    BigInt gcd = GCD( tmp, n );
    if( gcd > one && gcd < n )
    {
        if( ctrl.progress )
            Output("Found stage-1 factor of ",gcd);
        return gcd;
    }

    if( separateOdd )
    {
        const double twoLog = Log(2.);
        unsigned twoExponent = unsigned(nLog/twoLog);
        for( unsigned i=0; i<twoExponent; ++i )
        {
            // a = a*a (mod n)
            a *= a;
            a %= n;

            // gcd = GCD( a-1, n );
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

    return gcd;
}

// NOTE: Returns the GCD of stage 1 and overwrites a with a power of a
template<typename TSieve,typename TSieveSmall>
BigInt StageTwo
( const BigInt& n,
        BigInt& a,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
        TSieve previousBound,
        TSieve newBound,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    const BigInt& one = BigIntOne();
    sieve.SetLowerBound( previousBound+1 );

    BigInt diffPower, tmp, gcd;

    Int delayCounter=1;
    TSieve p, pLast=sieve.NextPrime();
    std::map<TSieve,BigInt> diffPowers;
    while( pLast <= newBound )
    {
        p = sieve.NextPrime();
        TSieve diff = p - pLast;
        auto search = diffPowers.find( diff );

        if( search == diffPowers.end() )
        {
            PowMod( a, diff, n, diffPower );
            diffPowers.insert( std::make_pair(diff,diffPower) );
        }

        a *= diffPowers[diff];
        a %= n;

        if( delayCounter >= ctrl.gcdDelay2 )
        {
            // gcd = GCD( a-1, n );
            tmp = a;
            tmp -= 1;
            GCD( tmp, n, gcd );
            if( gcd > one && gcd < n )
            {
                if( ctrl.progress )
                    Output("Found stage-2 factor of ",gcd);
                return gcd;
            }
            delayCounter = 0;
        }
        ++delayCounter;
    }

    // If the last iteration did not perform a GCD due to the delay
    if( ctrl.gcdDelay2 > 1 && delayCounter != 1 )
    {
        tmp = a;
        tmp -= 1;
        GCD( tmp, n, gcd );
        if( gcd > one && gcd < n )
        {
            if( ctrl.progress )
                Output("Found stage-2 factor of ",gcd);
        }
    }

    return gcd;
}

template<typename TSieve,typename TSieveSmall>
BigInt FindFactor
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    const BigInt& one = BigIntOne();
    const BigInt& two = BigIntTwo();

    TSieve smooth1 = ctrl.smooth1;
    TSieve smooth2 = ctrl.smooth2;
    TSieve smooth1Bound = Max( 1000000000ULL, 8ULL*smooth1 );
    TSieve smooth2Bound = Max( 10000000000ULL, 8ULL*smooth2 );

    // Ensure that we do not have a GCD of n appear too many times despite
    // separately checking the powers of two
    bool separateOdd=false;
    Int maxGCDFailures=10;
    Int numGCDFailures=0;

    BigInt gcd, tmp, diffPower;
    while( true )
    {
        // Uniformly select a in (Z/(n))* \ {1}
        BigInt a = SampleUniform( two, n );
        while( GCD( a, n ) != one )
        {
            a = SampleUniform( two, n );
        }

        gcd = StageOne( n, a, sieve, separateOdd, smooth1, ctrl );
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
        else if( gcd == one )
        {
            if( ctrl.progress )
                Output("Proceeding to stage-2 with a=",a);
        }
        else
            return gcd;

        // Run stage-2
        // -----------
        gcd = StageTwo( n, a, sieve, smooth1, smooth2, ctrl );
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
        else if( gcd == one )
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
        else
            return gcd;
    }
}

template<typename TSieve,typename TSieveSmall>
BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    DynamicSieve<TSieve,TSieveSmall> sieve;
    return FindFactor( n, sieve, ctrl );
}

} // namespace pollard_pm1

template<typename TSieve,typename TSieveSmall>
vector<BigInt> PollardPMinusOne
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl )
{
    vector<BigInt> factors;
    BigInt nRem = n;

    if( !ctrl.avoidTrialDiv )
    {
        // Start with trial division
        auto tinyFactors = TrialDivision( n, ctrl.trialDivLimit );
        for( auto tinyFactor : tinyFactors )
        {
            factors.push_back( tinyFactor );
            nRem /= tinyFactor;
            if( ctrl.progress )
                Output("Removed tiny factor of ",tinyFactor);
        }
    }
    if( nRem <= BigInt(1) )
        return factors;

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
        for( const auto& subfactor : subfactors )
            factors.push_back( subfactor );
        nRem /= factor;
    }
    PopIndent();
    sort( factors.begin(), factors.end() );
    return factors;
}

template<typename TSieve,typename TSieveSmall>
vector<BigInt> PollardPMinusOne
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
