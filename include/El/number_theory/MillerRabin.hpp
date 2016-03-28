/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_MILLER_RABIN_HPP
#define EL_NUMBER_THEORY_MILLER_RABIN_HPP

namespace El {

#ifdef EL_HAVE_MPC

// See Algorithm 8.2.2 from Henri Cohen's 
// "A course in computational algebraic number theory"
inline Primality MillerRabin( const BigInt& n, Int numReps )
{
    BigInt one(1), two(2);
    if( n <= one )
        LogicError("Miller-Rabin is invalid for n <= 1");
    else if( n == two )
        return PRIME;
    else if( Mod(n,two) == BigInt(0) )
        return COMPOSITE;

    // Decompose n-1 as 2^t*q, where q is odd
    // --------------------------------------
    BigInt nm1(n);
    nm1 -= 1;
    BigInt q;
    auto t =
      mpz_remove( q.Pointer(), nm1.LockedPointer(), two.LockedPointer() );

    BigInt a, b;
    for( Int c=0; c<numReps; ++c )
    {
        // Sample 1 < a < n-1, which is a subset of (Z/nZ)* if n is prime.
        // We leave out 1 since, obviously, 1^k = 1 for all k >= 0, and 
        // n-1 since (n-1)^2 = 1 (mod n), regardless of whether n is prime,
        // and so we would have (n-1)^k = +-1 (mod n) for all values of k and
        // the test would not be of any use.
        a = SampleUniform( two, nm1 );

        // b := a^q (mod n)
        PowMod( a, q, n, b );
        if( b == one )
        {
            // n is a strong probable prime to base a, as
            // a^(n-1) = (a^q)^(2^t) = 1^(2^t) = 1 (mod n)
            continue;
        }

        for( decltype(t) e=0; e<t-1; ++e )
        {
            if( b == nm1 )
            {
                // b^2 (mod n) will be equal to one, so a^(n-1) (mod n) = 1
                // and a is a strong probable prime to base a
                break;
            }
            b *= b;
            b %= n;

            if( b == one )
            {
                // The only square-roots of 1 (modulo a prime) are +-1, and
                // none of the square-roots that we have previously encountered
                // have been of this form, so n cannot be a prime.
                return COMPOSITE;
            }
        }
        if( b != nm1 )
        {
            // Since we did not break the above loop, b = a^(n-1)/2 (mod n),
            // and b != +-1 (mod n), so a^(n-1) != 1 (mod n) and n cannot be 
            // prime.
            return COMPOSITE;
        }
    }
    return PROBABLY_PRIME;
}

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_MILLER_RABIN_HPP
