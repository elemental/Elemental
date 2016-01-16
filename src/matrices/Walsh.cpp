/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T> 
void Walsh( Matrix<T>& A, Int k, bool binary )
{
    DEBUG_ONLY(CSE cse("Walsh"))
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, n );

    // Run a simple O(n^2 log n) algorithm for computing the entries
    // based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    auto walshFill = 
      [&]( Int i, Int j ) -> T
      {
        // Recurse on the quadtree, flipping the sign of the entry each
        // time we are in the bottom-right quadrant
        Unsigned r = (Unsigned)i;     
        Unsigned s = (Unsigned)j;
        Unsigned t = n;
        bool on = true;
        while( t != 1u )
        {
            t >>= 1;
            if( r >= t && s >= t )
                on = !on;
            r %= t;
            s %= t;
        }
        return ( on ? onValue : offValue );
      };
    IndexDependentFill( A, function<T(Int,Int)>(walshFill) );
}

template<typename T>
void Walsh( AbstractDistMatrix<T>& A, Int k, bool binary )
{
    DEBUG_ONLY(CSE cse("Walsh"))
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, n );

    // Run a simple O(n^2 log n) algorithm for computing the entries
    // based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    auto walshFill = 
      [&]( Int i, Int j ) -> T
      {
        // Recurse on the quadtree, flipping the sign of the entry each
        // time we are in the bottom-right quadrant
        Unsigned r = (Unsigned)i;     
        Unsigned s = (Unsigned)j;
        Unsigned t = n;
        bool on = true;
        while( t != 1u )
        {
            t >>= 1;
            if( r >= t && s >= t )
                on = !on;
            r %= t;
            s %= t;
        }
        return ( on ? onValue : offValue );
      };
    IndexDependentFill( A, function<T(Int,Int)>(walshFill) );
}

#define PROTO(T) \
  template void Walsh( Matrix<T>& A, Int k, bool binary ); \
  template void Walsh( AbstractDistMatrix<T>& A, Int k, bool binary );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
