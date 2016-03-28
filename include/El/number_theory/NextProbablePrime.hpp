/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_NEXT_PROBABLE_PRIME_HPP
#define EL_NUMBER_THEORY_NEXT_PROBABLE_PRIME_HPP

namespace El {

#ifdef EL_HAVE_MPC

// TODO: Use a custom implementation that incorporates 'numReps'
inline void NextProbablePrime( const BigInt& n, BigInt& nextPrime, Int numReps )
{
    mpz_nextprime( nextPrime.Pointer(), n.LockedPointer() );
}

inline BigInt NextProbablePrime( const BigInt& n, Int numReps )
{
    BigInt nextPrime;
    NextProbablePrime( n, nextPrime, numReps );
    return nextPrime;
}

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_NEXT_PROBABLE_PRIME_HPP
