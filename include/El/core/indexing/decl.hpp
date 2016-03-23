/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INDEXING_DECL_HPP
#define EL_INDEXING_DECL_HPP

namespace El {

// Indexing for element-wise distributions
// =======================================
Int Length( Int n, Int rank, Int firstRank, Int numProcs );
Int Length_( Int n, Int rank, Int firstRank, Int numProcs ) EL_NO_EXCEPT;

Int Length( Int n, Int shift, Int numProcs );
Int Length_( Int n, Int shift, Int numProcs ) EL_NO_EXCEPT;

Int MaxLength( Int n, Int numProcs );
Int MaxLength_( Int n, Int numProcs ) EL_NO_EXCEPT;

Int GlobalIndex( Int iLoc, Int shift, Int numProcs );

// Indexing for block distributions
// ================================
Int BlockedLength
( Int n, Int rank, Int firstRank, Int bsize, Int cut, Int numProcs );
Int BlockedLength_
( Int n, Int rank, Int firstRank, Int bsize, Int cut, Int numProcs )
EL_NO_EXCEPT;

Int BlockedLength ( Int n, Int shift, Int bsize, Int cut, Int numProcs );
Int BlockedLength_( Int n, Int shift, Int bsize, Int cut, Int numProcs )
EL_NO_EXCEPT;

Int MaxBlockedLength ( Int n, Int bsize, Int cut, Int numProcs );
Int MaxBlockedLength_( Int n, Int bsize, Int cut, Int numProcs ) EL_NO_EXCEPT;

Int GlobalBlockedIndex( Int iLoc, Int shift, Int bsize, Int cut, Int numProcs );

// Miscellaneous indexing routines
// ===============================

// Generalization of "%" operator which handles negative a in a way which
// still returns a result in [0,b). Note that b is assumed to be non-negative.
Int Mod( Int a, Int b );
Int Mod_( Int a, Int b ) EL_NO_EXCEPT;
#ifdef EL_HAVE_MPC
BigInt Mod( const BigInt& a, const BigInt& b );
BigInt Mod( const BigInt& a, const unsigned& b );
BigInt Mod( const BigInt& a, const unsigned long& b );
BigInt Mod_( const BigInt& a, const BigInt& b );
#endif

Int Shift( Int rank, Int firstRank, Int numProcs );
Int Shift_( Int rank, Int firstRank, Int numProcs ) EL_NO_EXCEPT;

Int LastOffset( Int n, Int bsize );

Int DiagonalLength( Int height, Int width, Int offset=0 ) EL_NO_EXCEPT;

Int GCD( Int a, Int b ); 
Int GCD_( Int a, Int b ) EL_NO_EXCEPT; 

#ifdef EL_HAVE_MPC
BigInt GCD( const BigInt& a, const BigInt& b );
#endif

Unsigned FlooredLog2( Unsigned n );
bool PowerOfTwo( Unsigned n );

#ifdef EL_HAVE_MPC
BigInt PowMod( const BigInt& base, const BigInt& exp, const BigInt& mod );
#endif

enum Primality
{
  PRIME,
  PROBABLY_PRIME,
  PROBABLY_COMPOSITE,
  COMPOSITE
};

#ifdef EL_HAVE_MPC
// Use a combination of trial divisions and Miller-Rabin 
// (with numReps representatives) to test for primality
Primality PrimalityTest( const BigInt& n, int numReps=20 );

// Return the first prime greater than n (with high likelihood)
BigInt NextPrime( const BigInt& n );
#endif

} // namespace El

#endif // ifndef EL_INDEXING_DECL_HPP
