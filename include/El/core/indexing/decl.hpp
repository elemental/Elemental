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
// TODO: ExtendedGCD
// TODO: LCM
// TODO: InvertMod

#ifdef EL_HAVE_MPC
BigInt GCD( const BigInt& a, const BigInt& b );
// A version which acts in-place to avoid an unnecessary allocation
void GCD( const BigInt& a, const BigInt& b, BigInt& gcd );

// Set gcd to be GCD(a,b) and (s,t) such that a*s + b*t = gcd.
void ExtendedGCD
( const BigInt& a, const BigInt& b, BigInt& gcd, BigInt& s, BigInt& t );

BigInt LCM( const BigInt& a, const BigInt& b );
void LCM( const BigInt& a, const BigInt& b, BigInt& lcm );

BigInt InvertMod( const BigInt& a, const BigInt& n );
void InvertMod( const BigInt& a, const BigInt& n, BigInt& aInv );
#endif

template<typename T,typename=EnableIf<IsIntegral<T>>>
void SqrtRem( const T& alpha, T& alphaSqrt, T& remainder );
#ifdef EL_HAVE_MPC
template<> void SqrtRem
( const BigInt& alpha, BigInt& alphaSqrt, BigInt& remainder );
#endif

template<typename T,typename=EnableIf<IsIntegral<T>>>
bool IsPerfectSquare( const T& alpha );
#ifdef EL_HAVE_MPC
template<> bool IsPerfectSquare( const BigInt& alpha );
#endif

Unsigned FlooredLog2( Unsigned n );
bool PowerOfTwo( Unsigned n );

#ifdef EL_HAVE_MPC
BigInt PowMod
( const BigInt& base,
  const BigInt& exp,
  const BigInt& mod );
void PowMod
( const BigInt& base,
  const BigInt& exp,
  const BigInt& mod,
        BigInt& result );

BigInt PowMod
( const BigInt& base,
        unsigned long exp,
  const BigInt& mod );
void PowMod
( const BigInt& base,
        unsigned long exp,
  const BigInt& mod,
        BigInt& result );
BigInt PowMod
( const BigInt& base,
        unsigned long long exp,
  const BigInt& mod );
void PowMod
( const BigInt& base,
        unsigned long long exp,
  const BigInt& mod,
        BigInt& result );
#endif

} // namespace El

#endif // ifndef EL_INDEXING_DECL_HPP
