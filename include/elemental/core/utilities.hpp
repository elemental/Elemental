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

int DiagonalLength( int height, int width, int offset=0 );

int GCD( int a, int b ); 

int RawGCD( int a, int b ); 

int LocalLength( int n, int shift, int modulus );

int RawLocalLength( int n, int shift, int modulus );

int LocalLength( int n, int index, int alignment, int modulus );

int RawLocalLength( int n, int index, int alignment, int modulus );

int MaxLocalLength( int n, int modulus );

int RawMaxLocalLength( int n, int modulus );

int Shift( int index, int alignment, int modulus );

int RawShift( int index, int alignment, int modulus );

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline int
elemental::DiagonalLength( int height, int width, int offset )
{
    if( offset > 0 )
    {
        const int remainingWidth = std::max(width-offset,0);
        return std::min(height,remainingWidth);
    }
    else
    {
        const int remainingHeight = std::max(height+offset,0);
        return std::min(remainingHeight,width);
    }
}

inline int
elemental::GCD( int a, int b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        throw std::logic_error( "GCD called with negative argument." );
#endif
    return elemental::RawGCD( a, b );
}

inline int
elemental::RawGCD( int a, int b )
{
    if( b == 0 )
        return a;
    else
        return RawGCD( b, a-b*(a/b) );
}

inline int
elemental::LocalLength( int n, int shift, int modulus )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
    if( n < 0 )
        throw std::logic_error( "n must be non-negative." );
    if( shift < 0 || shift >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", modulus=" << modulus;
        throw std::logic_error( msg.str() );
    }
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return elemental::RawLocalLength( n, shift, modulus );
}

inline int
elemental::RawLocalLength( int n, int shift, int modulus )
{
    return ( n > shift ? (n - shift - 1)/modulus + 1 : 0 );
}

inline int
elemental::LocalLength
( int n, int index, int alignment, int modulus )
{
#ifndef RELEASE
    PushCallStack("LocalLength");
#endif
    int shift = Shift( index, alignment, modulus );
    int localLength = LocalLength( n, shift, modulus );
#ifndef RELEASE
    PopCallStack();
#endif
    return localLength;
}

inline int
elemental::RawLocalLength
( int n, int index, int alignment, int modulus )
{
    int shift = RawShift( index, alignment, modulus );
    int localLength = RawLocalLength( n, shift, modulus );
    return localLength;
}

inline int
elemental::MaxLocalLength( int n, int modulus )
{
#ifndef RELEASE
    PushCallStack("MaxLocalLength");
    if( n < 0 )
        throw std::logic_error( "n must be non-negative." );
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return elemental::RawMaxLocalLength( n, modulus );
}

inline int
elemental::RawMaxLocalLength( int n, int modulus )
{
    return ( n > 0 ? (n - 1)/modulus + 1 : 0 );
}

// For determining the first global element of process row/column 
// 'index', with distribution alignment 'alignment' and number of process 
// rows/cols 'modulus'
inline int
elemental::Shift( int index, int alignment, int modulus )
{
#ifndef RELEASE
    PushCallStack("Shift");
    if( index < 0 || index >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid index: "
            << "index=" << index << ", modulus=" << modulus;
        throw std::logic_error( msg.str() );
    }
    if( alignment < 0 || alignment >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "alignment=" << alignment << ", modulus=" << modulus;
        throw std::logic_error( msg.str() );
    }
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return elemental::RawShift( index, alignment, modulus );
}

inline int
elemental::RawShift( int index, int alignment, int modulus )
{
    return (index + modulus - alignment) % modulus;
}

#endif /* ELEMENTAL_UTILITIES_HPP */

