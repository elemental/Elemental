/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_INDEXING_IMPL_HPP
#define EL_INDEXING_IMPL_HPP

namespace El {

// Indexing for element-wise distributions
// =======================================

inline Int Length( Int n, Int shift, Int stride )
{
    DEBUG_ONLY(
      CSE cse("Length");
      if( n < 0 )
          LogicError("n must be non-negative");
      if( shift < 0 || shift >= stride )
          LogicError("Invalid shift: shift=",shift,", stride=",stride);
      if( stride <= 0 )
          LogicError("Modulus must be positive");
    )
    return Length_( n, shift, stride );
}

inline Int Length_( Int n, Int shift, Int stride ) EL_NO_EXCEPT
{
    return ( n > shift ? (n - shift - 1)/stride + 1 : 0 );
}

inline Int
Length( Int n, Int rank, Int align, Int stride )
{
    DEBUG_ONLY(CSE cse("Length"))
    const Int shift = Shift( rank, align, stride );
    return Length( n, shift, stride );
}

inline Int Length_
( Int n, Int rank, Int align, Int stride ) EL_NO_EXCEPT
{
    const Int shift = Shift_( rank, align, stride );
    return Length_( n, shift, stride );
}

inline Int MaxLength( Int n, Int stride )
{
    DEBUG_ONLY(
      CSE cse("MaxLength");
      if( n < 0 )
          LogicError("n must be non-negative");
      if( stride <= 0 )
          LogicError("Modulus must be positive");
    )
    return MaxLength_( n, stride );
}

inline Int MaxLength_( Int n, Int stride ) EL_NO_EXCEPT
{ return ( n>0 ? (n-1)/stride + 1 : 0 ); }

inline Int GlobalIndex( Int iLoc, Int shift, Int numProcs )
{ return shift + iLoc*numProcs; }

// Indexing for block distributions
// ================================

inline Int BlockedLength( Int n, Int shift, Int bsize, Int cut, Int stride )
{
    DEBUG_ONLY(
      CSE cse("BlockedLength");
      if( n < 0 )
          LogicError("n must be non-negative");
      if( shift < 0 || shift >= stride )
          LogicError("Invalid shift: shift=",shift,", stride=",stride);
      if( stride <= 0 )
          LogicError("Modulus must be positive");
      // TODO: bsize and cut checks
    )
    return BlockedLength_( n, shift, bsize, cut, stride );
}

inline Int BlockedLength_
( Int n, Int shift, Int bsize, Int cut, Int stride ) EL_NO_EXCEPT
{
    Int length=0;

    // Handle the first block
    // ======================
    const Int firstLeftover = Min(n,bsize-cut);
    if( shift == 0 )
        length += firstLeftover;
    n -= firstLeftover;
    // Cycle each process's first block left one
    shift = Mod(shift-1,stride);

    // Handle the middle blocks
    // ========================
    const Int nBlock = n/bsize;
    const Int lengthBlock = Length_( nBlock, shift, stride );
    length += lengthBlock*bsize;
    n -= nBlock*bsize;
    // Cycle each process's first block left by nBlock
    shift = Mod(shift-nBlock,stride);

    // Handle the (possibly empty) last block
    // ======================================
    if( shift == 0 )
        length += n;

    return length;
}

inline Int
BlockedLength( Int n, Int rank, Int align, Int bsize, Int cut, Int stride )
{
    DEBUG_ONLY(CSE cse("BlockedLength"))
    const Int shift = Shift( rank, align, stride );
    return BlockedLength( n, shift, bsize, cut, stride );
}

inline Int BlockedLength_
( Int n, Int rank, Int align, Int bsize, Int cut, Int stride ) EL_NO_EXCEPT
{
    const Int shift = Shift_( rank, align, stride );
    return BlockedLength_( n, shift, bsize, cut, stride );
}

inline Int MaxBlockedLength( Int n, Int bsize, Int cut, Int stride )
{
    DEBUG_ONLY(
      CSE cse("MaxBlockedLength");
      if( n < 0 )
          LogicError("n must be non-negative");
      if( stride <= 0 )
          LogicError("Modulus must be positive");
    )
    return MaxBlockedLength_( n, bsize, cut, stride );
}

inline Int MaxBlockedLength_
( Int n, Int bsize, Int cut, Int stride ) EL_NO_EXCEPT
{ return BlockedLength_( n, 0, bsize, cut, stride ); }

inline Int 
GlobalBlockedIndex( Int iLoc, Int shift, Int bsize, Int cut, Int numProcs )
{ 
    // The number of global entries before the first block this process owns
    // data in begins (NOTE: this is negative if we own the first block and
    // the cut is nonzero)
    const Int iBefore = shift*bsize - cut;

    const Int iLocAdj = ( shift==0 ? iLoc+cut : iLoc );
    const Int numFilledLocalBlocks = iLocAdj / bsize;
    const Int iMid = numFilledLocalBlocks*bsize*numProcs;

    const Int iPost = iLocAdj-numFilledLocalBlocks*bsize;

    return iBefore + iMid + iPost;
}

// Miscellaneous indexing routines
// ===============================

inline Int Mod( Int a, Int b )
{
    DEBUG_ONLY(
      CSE cse("Mod");
      if( b <= 0 )
          LogicError("b is assumed to be positive");
    )
    return Mod_( a, b );
}

inline Int Mod_( Int a, Int b ) EL_NO_EXCEPT
{
    const Int rem = a % b;
    return ( rem >= 0 ? rem : rem+b );
}

// For determining the first index assigned to a given rank
inline Int Shift( Int rank, Int align, Int stride )
{
    DEBUG_ONLY(
      CSE cse("Shift");
      if( rank < 0 || rank >= stride )
          LogicError("Invalid rank: rank=",rank,", stride=",stride);
      if( align < 0 || align >= stride )
          LogicError("Invalid alignment: align=",align,", stride=",stride);
      if( stride <= 0 )
          LogicError("Stride must be positive");
    )
    return Shift_( rank, align, stride );
}

inline Int Shift_( Int rank, Int align, Int stride ) EL_NO_EXCEPT
{ return Mod(rank-align,stride); }


inline Int LastOffset( Int n, Int bsize )
{ return bsize*( Mod(n,bsize) ? n/bsize : (n/bsize)-1 ); }

inline Int
DiagonalLength( Int height, Int width, Int offset ) EL_NO_EXCEPT
{
    if( offset > 0 )
    {
        const Int remWidth = Max(width-offset,0);
        return Min(height,remWidth);
    }
    else
    {
        const Int remHeight = Max(height+offset,0);
        return Min(remHeight,width);
    }
}

inline Int GCD( Int a, Int b )
{
    DEBUG_ONLY(
      CSE cse("GCD");
      if( a < 0 || b < 0 )
          LogicError("GCD called with negative argument");
    )
    return GCD_( a, b );
}

inline Int GCD_( Int a, Int b ) EL_NO_EXCEPT
{
    if( b == 0 )
        return a;
    else
        return GCD_( b, a-b*(a/b) );
}

inline bool PowerOfTwo( Unsigned n )
{ return n && !(n & (n-1)); }

inline Unsigned Log2( Unsigned n )
{
    Unsigned result=0;
    while( n >>= 1 )
      ++result;
    return result;
}

} // namespace El

#endif // ifndef EL_INDEXING_IMPL_HPP
