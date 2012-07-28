/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
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
inline Int LocalLength( Int n, Int shift, Int numProcs )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( shift < 0 || shift >= numProcs )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", numProcs=" << numProcs;
        throw std::logic_error( msg.str().c_str() );
    }
    if( numProcs <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawLocalLength( n, shift, numProcs );
}

template<typename Int>
inline Int RawLocalLength( Int n, Int shift, Int numProcs )
{
    return ( n > shift ? (n - shift - 1)/numProcs + 1 : 0 );
}

template<typename Int>
inline Int 
LocalLength( Int n, Int rank, Int firstRank, Int numProcs )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
#endif
    Int shift = Shift( rank, firstRank, numProcs );
    Int localLength = LocalLength( n, shift, numProcs );
#ifndef RELEASE
    PopCallStack();
#endif
    return localLength;
}

template<typename Int>
inline Int RawLocalLength
( Int n, Int rank, Int firstRank, Int numProcs )
{
    Int shift = RawShift( rank, firstRank, numProcs );
    Int localLength = RawLocalLength( n, shift, numProcs );
    return localLength;
}

template<typename Int>
inline Int MaxLocalLength( Int n, Int numProcs )
{
#ifndef RELEASE
    PushCallStack("MaxLocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( numProcs <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawMaxLocalLength( n, numProcs );
}

template<typename Int>
inline Int RawMaxLocalLength( Int n, Int numProcs )
{
    return ( n > 0 ? (n - 1)/numProcs + 1 : 0 );
}

// For determining the first global element of process row/column 
// 'rank', with distribution alignment 'firstRank' and number of process 
// rows/cols 'numProcs'
template<typename Int>
inline Int Shift( Int rank, Int firstRank, Int numProcs )
{
#ifndef RELEASE
    PushCallStack("Shift");
    if( rank < 0 || rank >= numProcs )
    {
        std::ostringstream msg;
        msg << "Invalid rank: "
            << "rank=" << rank << ", numProcs=" << numProcs;
        throw std::logic_error( msg.str().c_str() );
    }
    if( firstRank < 0 || firstRank >= numProcs )
    {
        std::ostringstream msg;
        msg << "Invalid first rank: "
            << "first rank=" << firstRank << ", numProcs=" << numProcs;
        throw std::logic_error( msg.str().c_str() );
    }
    if( numProcs <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawShift( rank, firstRank, numProcs );
}

template<typename Int>
inline Int RawShift( Int rank, Int firstRank, Int numProcs )
{
    return (rank + numProcs - firstRank) % numProcs;
}

} // namespace elem
