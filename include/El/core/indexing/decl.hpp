/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
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

Int Shift( Int rank, Int firstRank, Int numProcs );
Int Shift_( Int rank, Int firstRank, Int numProcs ) EL_NO_EXCEPT;

Int LastOffset( Int n, Int bsize );

Int DiagonalLength( Int height, Int width, Int offset=0 ) EL_NO_EXCEPT;

Int GCD( Int a, Int b ); 
Int GCD_( Int a, Int b ) EL_NO_EXCEPT; 

Unsigned Log2( Unsigned n );
bool PowerOfTwo( Unsigned n );

} // namespace El

#endif // ifndef EL_INDEXING_DECL_HPP
