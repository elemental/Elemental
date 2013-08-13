/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_INDEXING_IMPL_HPP
#define ELEM_CORE_INDEXING_IMPL_HPP

namespace elem {

inline Int
DiagonalLength( Int height, Int width, Int offset )
{
    if( offset > 0 )
    {
        Int remainingWidth = Max(width-offset,0);
        return Min(height,remainingWidth);
    }
    else
    {
        Int remainingHeight = Max(height+offset,0);
        return Min(remainingHeight,width);
    }
}

inline Int GCD( Int a, Int b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        LogicError("GCD called with negative argument");
#endif
    return GCD_( a, b );
}

inline Int GCD_( Int a, Int b )
{
    if( b == 0 )
        return a;
    else
        return GCD_( b, a-b*(a/b) );
}

inline Int Length( Int n, Int shift, Int stride )
{
#ifndef RELEASE
    CallStackEntry entry("Length");
    if( n < 0 )
        LogicError("n must be non-negative");
    if( shift < 0 || shift >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", stride=" << stride;
        LogicError( msg.str() );
    }
    if( stride <= 0 )
        LogicError("Modulus must be positive");
#endif
    return Length_( n, shift, stride );
}

inline Int Length_( Int n, Int shift, Int stride )
{
    return ( n > shift ? (n - shift - 1)/stride + 1 : 0 );
}

inline Int
Length( Int n, Int rank, Int alignment, Int stride )
{
#ifndef RELEASE
    CallStackEntry entry("Length");
#endif
    Int shift = Shift( rank, alignment, stride );
    return Length( n, shift, stride );
}

inline Int Length_
( Int n, Int rank, Int alignment, Int stride )
{
    Int shift = Shift_( rank, alignment, stride );
    return Length_( n, shift, stride );
}

inline Int MaxLength( Int n, Int stride )
{
#ifndef RELEASE
    CallStackEntry entry("MaxLength");
    if( n < 0 )
        LogicError("n must be non-negative");
    if( stride <= 0 )
        LogicError("Modulus must be positive");
#endif
    return MaxLength_( n, stride );
}

inline Int MaxLength_( Int n, Int stride )
{
    return ( n > 0 ? (n - 1)/stride + 1 : 0 );
}

// For determining the first index assigned to a given rank
inline Int Shift( Int rank, Int alignment, Int stride )
{
#ifndef RELEASE
    CallStackEntry entry("Shift");
    if( rank < 0 || rank >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid rank: "
            << "rank=" << rank << ", stride=" << stride;
        LogicError( msg.str() );
    }
    if( alignment < 0 || alignment >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "alignment=" << alignment << ", stride=" << stride;
        LogicError( msg.str() );
    }
    if( stride <= 0 )
        LogicError("Stride must be positive");
#endif
    return Shift_( rank, alignment, stride );
}

inline Int Shift_( Int rank, Int alignment, Int stride )
{ return (rank + stride - alignment) % stride; }

} // namespace elem

#endif // ifndef ELEM_CORE_INDEXING_IMPL_HPP
