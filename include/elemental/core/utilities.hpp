/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_UTILITIES_HPP
#define ELEMENTAL_UTILITIES_HPP 1

namespace elemental {

template<typename Integer>
Integer DiagonalLength( Integer height, Integer width, Integer offset=0 );

template<typename Integer>
Integer GCD( Integer a, Integer b ); 

template<typename Integer>
Integer RawGCD( Integer a, Integer b ); 

template<typename Integer>
Integer LocalLength( Integer n, Integer shift, Integer modulus );

template<typename Integer>
Integer RawLocalLength( Integer n, Integer shift, Integer modulus );

template<typename Integer>
Integer LocalLength
( Integer n, Integer index, Integer alignment, Integer modulus );

template<typename Integer>
Integer RawLocalLength
( Integer n, Integer index, Integer alignment, Integer modulus );

template<typename Integer>
Integer MaxLocalLength( Integer n, Integer modulus );

template<typename Integer>
Integer RawMaxLocalLength( Integer n, Integer modulus );

template<typename Integer>
Integer Shift( Integer index, Integer alignment, Integer modulus );

template<typename Integer>
Integer RawShift( Integer index, Integer alignment, Integer modulus );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename Integer>
inline Integer 
DiagonalLength( Integer height, Integer width, Integer offset )
{
    if( offset > 0 )
    {
        Integer remainingWidth = std::max(width-offset,0);
        return std::min(height,remainingWidth);
    }
    else
    {
        Integer remainingHeight = std::max(height+offset,0);
        return std::min(remainingHeight,width);
    }
}

template<typename Integer>
inline Integer GCD( Integer a, Integer b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        throw std::logic_error("GCD called with negative argument");
#endif
    return RawGCD( a, b );
}

template<typename Integer>
inline Integer RawGCD( Integer a, Integer b )
{
    if( b == 0 )
        return a;
    else
        return RawGCD( b, a-b*(a/b) );
}

template<typename Integer>
inline Integer LocalLength( Integer n, Integer shift, Integer modulus )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( shift < 0 || shift >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", modulus=" << modulus;
        throw std::logic_error( msg.str().c_str() );
    }
    if( modulus <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawLocalLength( n, shift, modulus );
}

template<typename Integer>
inline Integer RawLocalLength( Integer n, Integer shift, Integer modulus )
{
    return ( n > shift ? (n - shift - 1)/modulus + 1 : 0 );
}

template<typename Integer>
inline Integer 
LocalLength( Integer n, Integer index, Integer alignment, Integer modulus )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
#endif
    Integer shift = Shift( index, alignment, modulus );
    Integer localLength = LocalLength( n, shift, modulus );
#ifndef RELEASE
    PopCallStack();
#endif
    return localLength;
}

template<typename Integer>
inline Integer RawLocalLength
( Integer n, Integer index, Integer alignment, Integer modulus )
{
    Integer shift = RawShift( index, alignment, modulus );
    Integer localLength = RawLocalLength( n, shift, modulus );
    return localLength;
}

template<typename Integer>
inline Integer MaxLocalLength( Integer n, Integer modulus )
{
#ifndef RELEASE
    PushCallStack("MaxLocalLength");
    if( n < 0 )
        throw std::logic_error("n must be non-negative");
    if( modulus <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawMaxLocalLength( n, modulus );
}

template<typename Integer>
inline Integer RawMaxLocalLength( Integer n, Integer modulus )
{
    return ( n > 0 ? (n - 1)/modulus + 1 : 0 );
}

// For determining the first global element of process row/column 
// 'index', with distribution alignment 'alignment' and number of process 
// rows/cols 'modulus'
template<typename Integer>
inline Integer Shift( Integer index, Integer alignment, Integer modulus )
{
#ifndef RELEASE
    PushCallStack("Shift");
    if( index < 0 || index >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid index: "
            << "index=" << index << ", modulus=" << modulus;
        throw std::logic_error( msg.str().c_str() );
    }
    if( alignment < 0 || alignment >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "alignment=" << alignment << ", modulus=" << modulus;
        throw std::logic_error( msg.str().c_str() );
    }
    if( modulus <= 0 )
        throw std::logic_error("Modulus must be positive");
    PopCallStack();
#endif
    return RawShift( index, alignment, modulus );
}

template<typename Integer>
inline Integer RawShift( Integer index, Integer alignment, Integer modulus )
{
    return (index + modulus - alignment) % modulus;
}

} // namespace elemental

#endif /* ELEMENTAL_UTILITIES_HPP */

