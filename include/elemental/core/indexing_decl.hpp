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

template<typename I>
I DiagonalLength( I height, I width, I offset=0 );

template<typename I>
I GCD( I a, I b ); 
template<typename I>
I GCD_( I a, I b ); 

Int Length( Int n, Int shift, Int numProcs );
Int Length_( Int n, Int shift, Int numProcs );

Int Length( Int n, Int rank, Int firstRank, Int numProcs );
Int Length_( Int n, Int rank, Int firstRank, Int numProcs );

Int MaxLength( Int n, Int numProcs );
Int MaxLength_( Int n, Int numProcs );

Int Shift( Int rank, Int firstRank, Int numProcs );
Int Shift_( Int rank, Int firstRank, Int numProcs );

} // namespace elem

#endif // ifndef ELEM_CORE_INDEXING_DECL_HPP
