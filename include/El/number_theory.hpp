/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_HPP
#define EL_NUMBER_THEORY_HPP

namespace El {

template<typename T=unsigned long long,
         typename TSmall=unsigned>
struct DynamicSieve
{
public:
    DynamicSieve
    ( T lowerBound=3,
      bool keepAll=false,
      TSmall segmentSize=32768 );

    void SetLowerBound( T lowerBound );
    void SetStorage( bool keepAll );
    void Generate( T upperBound );

    T NextPrime();

    // We could use TSmall if keepAll was false
    vector<T> oddPrimes;

private:
    bool keepAll_;
    T lowerBound_; // always return the first prime >= lowerBound_
    T oddPrimeBound_;
    vector<TSmall> oddPrimeStarts_;

    // The segment offset will be odd (>=3) and the associated numbers are
    // the offset plus twice the index.
    vector<char> segmentTable_;
    TSmall segmentSize_;
    T segmentOffset_;
    TSmall segmentIndex_;
    // Attempt to update 'halvedIndex' until oddOffset + 2*halvedIndex is 
    // corresponds to a precomputed prime (i.e., table[halvedIndex] is one)
    bool SeekSegmentPrime();
    TSmall ComputeSegmentStart( const TSmall& p ) const;

    void AugmentPrimes( T numPrimes );
    void MoveSegmentOffset( T segmentOffset );

    T PrimeCountingEstimate( T n ) const;

    void SieveSegment();
    void FormNewSegment();
};

// For retrieving a global-scope sieve for trial division
DynamicSieve<unsigned long long,unsigned>& TrialDivisionSieve();

// Return the prime factors (up to the specified limit) found through trial div
vector<unsigned long long>
TrialDivision( unsigned long long n, unsigned long long limit=53 );
#ifdef EL_HAVE_MPC
vector<unsigned long long>
TrialDivision( const BigInt& n, unsigned long long limit=53 );
#endif

// Simply return whether or not a factor was found through trial division
bool HasTinyFactor( unsigned long long n, unsigned long long limit=53 );
#ifdef EL_HAVE_MPC
bool HasTinyFactor( const BigInt& n, unsigned long long limit=53 );
#endif

#ifdef EL_HAVE_MPC

unsigned long PowerDecomp
( const BigInt& n,
        BigInt& q,
  const BigInt& factor=BigIntTwo() );

BigInt SqrtModPrime( const BigInt& n, const BigInt& p );
void SqrtModPrime( const BigInt& n, const BigInt& p, BigInt& nSqrt );

int LegendreSymbol( const BigInt& n, const BigInt& p );
int JacobiSymbol( const BigInt& m, const BigInt& n );

enum Primality
{
  PRIME,
  PROBABLY_PRIME,
  PROBABLY_COMPOSITE,
  COMPOSITE
};

Primality MillerRabin( const BigInt& n, const BigInt& a=BigIntTwo() );
Primality MillerRabinSequence( const BigInt& n, Int numReps=30 );

// Use a combination of trial divisions and Miller-Rabin 
// (with numReps representatives) to test for primality.
Primality PrimalityTest( const BigInt& n, Int numReps=30 );

// Return the first prime greater than n (with high likelihood)
BigInt NextProbablePrime( const BigInt& n, Int numReps=30 );
void NextProbablePrime( const BigInt& n, BigInt& nextPrime, Int numReps=30 );

namespace factor {

struct PollardRhoCtrl
{
    Int a0=1;
    Int a1=-1;
    unsigned long numSteps=1u;
    BigInt x0=BigIntTwo();
    Int gcdDelay=100;

    // For trial division
    bool avoidTrialDiv=false;
    unsigned long long trialDivLimit=53ULL;

    // For Miller-Rabin primality testing
    Int numReps=30;

    bool progress=false;
    bool time=false;
};

vector<BigInt> PollardRho
( const BigInt& n,
  const PollardRhoCtrl& ctrl=PollardRhoCtrl() );

namespace pollard_rho {

BigInt FindDivisor
( const BigInt& n,
        Int a=1,
  const PollardRhoCtrl& ctrl=PollardRhoCtrl() );

} // namespace pollard_rho

template<typename TSieve=unsigned long long>
struct PollardPMinusOneCtrl
{
    // Stage one
    TSieve smooth1=TSieve(1000000ULL);
    bool jumpstart1=false;
    TSieve start1=2;

    // Stage two
    TSieve smooth2=TSieve(10000000ULL);
    Int gcdDelay2=100;

    // For trial division
    bool avoidTrialDiv=false;
    unsigned long long trialDivLimit=53ULL;

    // For Miller-Rabin primality testing
    Int numReps=30;

    bool progress=false;
    bool time=false;

    bool checkpoint=false;
    Int checkpointFreq=1000000;
};

template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl<TSieve>& ctrl=
        PollardPMinusOneCtrl<TSieve>() );

template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
vector<BigInt> PollardPMinusOne
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl=
        PollardPMinusOneCtrl<TSieve>() );

namespace pollard_pm1 {

template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl<TSieve>& ctrl=
        PollardPMinusOneCtrl<TSieve>() );

template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
BigInt FindFactor
( const BigInt& n,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
  const PollardPMinusOneCtrl<TSieve>& ctrl=
        PollardPMinusOneCtrl<TSieve>() );

template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
BigInt StageOne
( const BigInt& n,
        BigInt& a,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
        bool separateOdd,
        TSieve primeBound,
  const PollardPMinusOneCtrl<TSieve>& ctrl );
template<typename TSieve=unsigned long long,
         typename TSieveSmall=unsigned>
BigInt StageTwo
( const BigInt& n,
        BigInt& a,
        DynamicSieve<TSieve,TSieveSmall>& sieve,
        TSieve previousBound,
        TSieve newBound,
  const PollardPMinusOneCtrl<TSieve>& ctrl );

} // namespace pollard_pm1

} // namespace factor

bool IsPrimitiveRoot
( const BigInt& primitive,
  const BigInt& p,
  const vector<BigInt>& pm1Factors,
        bool progress=false );
bool IsPrimitiveRoot
( const BigInt& primitive,
  const BigInt& p,
        bool progress=false,
  const factor::PollardRhoCtrl& ctrl=factor::PollardRhoCtrl() );

// Return a primitive root of a prime number p
BigInt PrimitiveRoot( const BigInt& p, Int numReps=30 );
void PrimitiveRoot( const BigInt& p, BigInt& primitive, Int numReps=30 );

namespace dlog {

struct PollardRhoCtrl
{
    BigInt a0=0;
    BigInt b0=0;
    bool multistage=true;
    bool assumePrime=false;
    factor::PollardRhoCtrl factorCtrl;

    bool progress=false;
    bool time=false;
};

// Return k such that r^k = q (mod p)
BigInt PollardRho
( const BigInt& q,
  const BigInt& r,
  const BigInt& p,
  const PollardRhoCtrl& ctrl=PollardRhoCtrl() );

} // namespace dlog

#endif // ifdef EL_HAVE_MPC

} // namespace El

#include <El/number_theory/DynamicSieve.hpp>
#include <El/number_theory/TrialDivision.hpp>

#include <El/number_theory/PowerDecomp.hpp>
#include <El/number_theory/SqrtModPrime.hpp>
#include <El/number_theory/LegendreSymbol.hpp>
#include <El/number_theory/JacobiSymbol.hpp>
#include <El/number_theory/MillerRabin.hpp>
#include <El/number_theory/PrimalityTest.hpp>
#include <El/number_theory/NextProbablePrime.hpp>
#include <El/number_theory/factor/PollardRho.hpp>
#include <El/number_theory/factor/PollardPMinusOne.hpp>
#include <El/number_theory/PrimitiveRoot.hpp>
#include <El/number_theory/dlog/PollardRho.hpp>

#include <El/number_theory/lattice.hpp>

#endif // ifndef EL_NUMBER_THEORY_HPP
