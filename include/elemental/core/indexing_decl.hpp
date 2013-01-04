/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename Int>
Int DiagonalLength( Int height, Int width, Int offset=0 );

template<typename Int>
Int GCD( Int a, Int b ); 

template<typename Int>
Int RawGCD( Int a, Int b ); 

template<typename Int>
Int LocalLength( Int n, Int shift, Int numProcs );

template<typename Int>
Int RawLocalLength( Int n, Int shift, Int numProcs );

template<typename Int>
Int LocalLength
( Int n, Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int RawLocalLength
( Int n, Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int MaxLocalLength( Int n, Int numProcs );

template<typename Int>
Int RawMaxLocalLength( Int n, Int numProcs );

template<typename Int>
Int Shift( Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int RawShift( Int rank, Int firstRank, Int numProcs );

} // namespace elem
