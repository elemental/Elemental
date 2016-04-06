/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_PRIMITIVE_ROOT_HPP
#define EL_NUMBER_THEORY_PRIMITIVE_ROOT_HPP

namespace El {

#ifdef EL_HAVE_MPC

inline bool IsPrimitiveRoot
( const BigInt& primitive,
  const BigInt& p,
  const vector<BigInt>& pm1Factors,
        bool progress )
{
    BigInt one(1), two(2);
    BigInt primMod(primitive);
    primMod %= p;

    // Try to early exit
    if( p == two )
        return primMod == one;

    if( IsPerfectSquare(primMod) )
        return false;

    Int numReps=30; // TODO: Make configurable
    Primality primality = PrimalityTest( p, numReps );
    if( primality == COMPOSITE )
        LogicError(p," was composite according to Miller-Rabin");

    BigInt pm1(p);
    pm1 -= 1;

    const Int numFactors = pm1Factors.size();
    BigInt exponent, primModPower;
    for( Int i=0; i<numFactors; ++i )
    {
        if( i > 0 && pm1Factors[i] == pm1Factors[i-1] )
            continue;
        exponent = pm1;
        exponent /= pm1Factors[i];
        PowMod( primMod, exponent, p, primModPower ); 
        if( primModPower == one )
        {
            if( progress )
                Output
                ("Primitive test failed for ",
                 primMod,"^",exponent," = ",primModPower," (mod ",p,")");
            return false;
        }
    }
    return true;
}

inline bool IsPrimitiveRoot
( const BigInt& primitive,
  const BigInt& p,
        bool progress,
  const factor::PollardRhoCtrl& ctrl )
{
    BigInt primMod(primitive);
    primMod %= p;

    // Try to early exit
    if( p == BigInt(2) )
        return primMod == BigInt(1);

    if( IsPerfectSquare(primMod) )
        return false;

    BigInt pm1(p);
    pm1 -= 1;
    auto pm1Factors = factor::PollardRho( pm1, ctrl );
    return IsPrimitiveRoot( primMod, p, pm1Factors, progress );
}

inline void PrimitiveRoot
( const BigInt& p,
        BigInt& primitive,
        int numReps,
  const factor::PollardRhoCtrl& ctrl )
{
    if( p == BigInt(2) )
    {
        primitive = 1;
        return;
    }

    // Test the (probable) primality with Miller-Rabin as a sanity check
    Primality primality = PrimalityTest( p, numReps );
    if( primality == COMPOSITE )
        LogicError(p," was composite according to Miller-Rabin");

    // Factor p-1 into (probable) primes using Pollard's rho
    BigInt pm1(p);
    pm1 -= 1;
    auto factors = factor::PollardRho( pm1, ctrl );

    // Find an a in Zp^* such that powers of a generate Zp^*
    // (Equivalently, check that a^((p-1)/p_i) = -1 for each unique prime 
    //  factor p_i of p-1)
    // TODO: Consider inlining IsPrimitiveRoot to avoid temporary allocations
    for( BigInt a=2; a<=pm1; ++a )
    {
        if( IsPrimitiveRoot( a, p, factors ) )
        {
            primitive = a;
            return;
        }
    }
    RuntimeError("Could not find a primitive of ",p);
}

inline BigInt PrimitiveRoot
( const BigInt& p,
        int numReps,
  const factor::PollardRhoCtrl& ctrl=factor::PollardRhoCtrl() )
{
    BigInt primitive;
    PrimitiveRoot( p, primitive, numReps, ctrl );
    return primitive;
}
    
#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_PRIMITIVE_ROOT_HPP
