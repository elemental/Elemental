/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_SQRT_MOD_PRIME_HPP
#define EL_NUMBER_THEORY_SQRT_MOD_PRIME_HPP

namespace El {

#ifdef EL_HAVE_MPC

// This is a simple implementation of Tonelli-Shanks as given in Algorithm
// 1.5.1 (Square Root Mod p) in Henri Cohen's
// "A course in computational algebraic number theory"
inline void SqrtModPrime( const BigInt& n, const BigInt& p, BigInt& x )
{
    BigInt one(1), two(2);
    if( p == two )
    {
        // Squaring is the identity operation in Z/2Z
        x = n;
        x %= two;
        return;
    }
    // TODO: Optionally ensure that p is prime?

    // Decompose p-1 as 2^e*q, where q is odd
    // --------------------------------------
    BigInt pm1(p);
    pm1 -= 1;
    BigInt q;
    auto e =
      mpz_remove( q.Pointer(), pm1.LockedPointer(), two.LockedPointer() ); 

    // Find a quadratic non-residue, a, and set z := a^q (mod p)
    // ---------------------------------------------------------
    BigInt a = SampleUniform( one, p );
    while( LegendreSymbol(a,p) != -1 )
    {
        a = SampleUniform( one, p );
    }
    BigInt z = PowMod( a, q, p ); 

    // Initialize
    // ----------
    auto r = e;
    BigInt y(z);
    // x := n^((q-1)/2) (mod p)
    x = PowMod( n, (q-1)/2, p );
    // b := n*x^2 (mod p)
    BigInt b(n);
    b *= x;
    b *= x;
    b %= p;
    // x := n*x (mod p)
    x *= n;
    x %= p;

    BigInt bPow, yExp, t;
    while( true )
    {
        // Find exponent
        // -------------
        if( b == one )
            return;
        bPow = b;
        bPow *= b;
        bPow %= p;
        decltype(r) m=1;
        for( ; m<r; ++m )
        {
            if( bPow == one )
                break;
            // NOTE: This is not needed if m==r-1
            bPow *= bPow;
            bPow %= p;
        }
        if( m == r )
            LogicError(n," is not a quadratic residue mod ",p);
    
        // Reduce exponent
        // ---------------
        Pow( two, r-m-1, yExp );
        PowMod( y, yExp, p, t );
        y = t;
        y *= t;
        y %= p;
        r = m;
        x *= t;
        x %= p;
        b *= y;
        b %= p;
    }
}

inline BigInt SqrtModPrime( const BigInt& n, const BigInt& p )
{
    BigInt nSqrt;
    SqrtModPrime( n, p, nSqrt );
    return nSqrt;
}

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_SQRT_MOD_PRIME_HPP
