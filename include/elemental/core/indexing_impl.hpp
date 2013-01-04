/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename Int>
inline Int 
DiagonalLength( Int height, Int width, Int offset )
{
    if( offset > 0 )
    {
        Int remainingWidth = std::max(width-offset,0);
        return std::min(height,remainingWidth);
    }
    else
    {
        Int remainingHeight = std::max(height+offset,0);
        return std::min(remainingHeight,width);
    }
}

template<typename Int>
inline Int GCD( Int a, Int b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        throw std::logic_error("GCD called with negative argument");
#endif
    return RawGCD( a, b );
}

template<typename Int>
inline Int RawGCD( Int a, Int b )
{
    if( b == 0 )
        return a;
    else
        return RawGCD( b, a-b*(a/b) );
}

template<typename Int>
inline Int LocalLength( Int n, Int shift, Int stride )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( shift < 0 || shift >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", stride=" << stride;
        throw std::logic_error( msg.str().c_str() );
    }
    if( stride <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawLocalLength( n, shift, stride );
}

template<typename Int>
inline Int RawLocalLength( Int n, Int shift, Int stride )
{
    return ( n > shift ? (n - shift - 1)/stride + 1 : 0 );
}

template<typename Int>
inline Int 
LocalLength( Int n, Int rank, Int alignment, Int stride )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
#endif
    Int shift = Shift( rank, alignment, stride );
    Int localLength = LocalLength( n, shift, stride );
#ifndef RELEASE
    PopCallStack();
#endif
    return localLength;
}

template<typename Int>
inline Int RawLocalLength
( Int n, Int rank, Int alignment, Int stride )
{
    Int shift = RawShift( rank, alignment, stride );
    Int localLength = RawLocalLength( n, shift, stride );
    return localLength;
}

template<typename Int>
inline Int MaxLocalLength( Int n, Int stride )
{
#ifndef RELEASE
    PushCallStack("MaxLocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( stride <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawMaxLocalLength( n, stride );
}

template<typename Int>
inline Int RawMaxLocalLength( Int n, Int stride )
{
    return ( n > 0 ? (n - 1)/stride + 1 : 0 );
}

// For determining the first index assigned to a given rank
template<typename Int>
inline Int Shift( Int rank, Int alignment, Int stride )
{
#ifndef RELEASE
    PushCallStack("Shift");
    if( rank < 0 || rank >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid rank: "
            << "rank=" << rank << ", stride=" << stride;
        throw std::logic_error( msg.str().c_str() );
    }
    if( alignment < 0 || alignment >= stride )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "alignment=" << alignment << ", stride=" << stride;
        throw std::logic_error( msg.str().c_str() );
    }
    if( stride <= 0 )
        throw std::logic_error("Stride must be positive");
    PopCallStack();
#endif
    return RawShift( rank, alignment, stride );
}

template<typename Int>
inline Int RawShift( Int rank, Int alignment, Int stride )
{
    return (rank + stride - alignment) % stride;
}

} // namespace elem
