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

// A dynamic version of the sieve of Eratosthenes
template<typename TUnsigned=unsigned long long>
struct DynamicSieve
{
public:
    DynamicSieve( TUnsigned lowerBound=3 );

    TUnsigned NextPrime();

    // A set of small odd primes for sieving larger primes
    vector<TUnsigned> oddPrimes;

private:
    vector<byte> table;
    TUnsigned oddOffset, halvedIndex;
    bool freshTable;

    // The smallest indices into the above table of each of the 'oddPrimes'
    // that could lead to said value of 'oddPrimes' being the smallest prime
    // factor
    vector<TUnsigned> starts;

    // Attempt to update 'halvedIndex' until oddOffset + 2*halvedIndex is 
    // corresponds to a precomputed prime (i.e., table[halvedIndex] is one)
    bool SeekPrecomputedPrime();

    TUnsigned ComputeStart( const TUnsigned& p ) const;

    void AugmentPrimes();
    void Sieve();
    void FormNewTable();

    TUnsigned CurrentPrime() const;
};

#ifdef EL_HAVE_MPC

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

Primality MillerRabin( const BigInt& n, Int numReps=30 );

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
    BigInt x0=BigInt(2);
    Int gcdDelay=100;
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

struct PollardPMinusOneCtrl
{
    unsigned long long smooth1=1000000ULL;
    unsigned long long smooth2=10000000ULL;
    Int numReps=30;
    bool progress=false;
    bool time=false;
};

vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() );

vector<BigInt> PollardPMinusOne
( const BigInt& n,
        DynamicSieve<unsigned long long>& sieve,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() );

namespace pollard_pm1 {

BigInt FindFactor
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() );

BigInt FindFactor
( const BigInt& n,
        DynamicSieve<unsigned long long>& sieve,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() );

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

#endif // ifndef EL_NUMBER_THEORY_HPP
