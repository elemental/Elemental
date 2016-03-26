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

#ifdef EL_HAVE_MPC

namespace factor {

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
    BigInt smoothness=BigInt(100000);
    int numReps=30;
    bool progress=false;
    bool time=false;
};

vector<BigInt> PollardPMinusOne
( const BigInt& n,
  const PollardPMinusOneCtrl& ctrl=PollardPMinusOneCtrl() );

namespace pollard_pm1 {

BigInt FindFactor
( const BigInt& n,
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
BigInt PrimitiveRoot( const BigInt& p, int numReps=30 );
void PrimitiveRoot( const BigInt& p, BigInt& primitive, int numReps=30 );

namespace dlog {

struct PollardRhoCtrl
{
    BigInt a0=0;
    BigInt b0=0;
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

#include <El/number_theory/factor/PollardRho.hpp>
#include <El/number_theory/factor/PollardPMinusOne.hpp>
#include <El/number_theory/PrimitiveRoot.hpp>
#include <El/number_theory/dlog/PollardRho.hpp>

#endif // ifndef EL_NUMBER_THEORY_HPP
