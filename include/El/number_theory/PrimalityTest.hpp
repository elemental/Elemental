/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_NUMBER_THEORY_PRIMALITY_TEST_HPP
#define EL_NUMBER_THEORY_PRIMALITY_TEST_HPP

namespace El {

#ifdef EL_HAVE_MPC

// TODO: A custom algorithm wrapping our Miller-Rabin
inline Primality PrimalityTest( const BigInt& n, Int numReps )
{
    int result = mpz_probab_prime_p( n.LockedPointer(), int(numReps) );
    if( result == 2 )
        return PRIME;
    else if( result == 1 )
        return PROBABLY_PRIME;
    else
        return COMPOSITE;
}

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_NUMBER_THEORY_PRIMALITY_TEST_HPP
