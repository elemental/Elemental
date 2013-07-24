/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_INDEXING_DECL_HPP
#define ELEM_CORE_INDEXING_DECL_HPP

namespace elem {

template<typename Int>
Int DiagonalLength( Int height, Int width, Int offset=0 );

template<typename Int>
Int GCD( Int a, Int b ); 
template<typename Int>
Int GCD_( Int a, Int b ); 

template<typename Int>
Int Length( Int n, Int shift, Int numProcs );
template<typename Int>
Int Length_( Int n, Int shift, Int numProcs );

template<typename Int>
Int Length( Int n, Int rank, Int firstRank, Int numProcs );
template<typename Int>
Int Length_( Int n, Int rank, Int firstRank, Int numProcs );

template<typename Int>
Int MaxLength( Int n, Int numProcs );
template<typename Int>
Int MaxLength_( Int n, Int numProcs );

template<typename Int>
Int Shift( Int rank, Int firstRank, Int numProcs );
template<typename Int>
Int Shift_( Int rank, Int firstRank, Int numProcs );

} // namespace elem

#endif // ifndef ELEM_CORE_INDEXING_DECL_HPP
