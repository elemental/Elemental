/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_UTILITIES_HPP
#define ELEMENTAL_UTILITIES_HPP 1

namespace Elemental
{
    namespace utilities
    {
        int
        GCD
        ( const int a, const int b ); 

        int
        LocalLength
        ( const int n, const int shift, const int modulus );

        int
        MaxLocalLength
        ( const int n, const int modulus );

        int
        Shift
        ( const int index, const int align, const int modulus );
    }
}

/*------------------------------------------------------------------------------
 * Utility function inlines                                
 *----------------------------------------------------------------------------*/

inline int
Elemental::utilities::GCD
( const int a, const int b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        throw "GCD called with negative argument.";
#endif
    if( b == 0 )
        return a;
    else
        return GCD( b, a-b*(a/b) );
}

inline int
Elemental::utilities::LocalLength
( const int n, const int shift, const int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::LocalLength");
    if( n < 0 )
        throw "n must be non-negative.";
    if( shift < 0 || shift >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", modulus=" << modulus << std::endl;
        const std::string s = msg.str();
        throw s.c_str();
    }
    if( modulus <= 0 )
        throw "Modulus must be positive.";
    PopCallStack();
#endif
    return ( n > shift ? (n - shift - 1)/modulus + 1 : 0 );
}

inline int
Elemental::utilities::MaxLocalLength
( const int n, const int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::MaxLocalLength");
    if( n < 0 )
        throw "n must be non-negative.";
    if( modulus <= 0 )
        throw "Modulus must be positive.";
    PopCallStack();
#endif
    return ( n > 0 ? (n - 1)/modulus + 1 : 0 );
}

// For determining the first global element of process row/column 
// 'index', with distribution alignment 'align' and number of process 
// rows/cols 'modulus'
inline int
Elemental::utilities::Shift
( const int index, const int align, const int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::Shift");
    if( index < 0 || index >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid index: "
            << "index=" << index << ", modulus=" << modulus << std::endl;
        const std::string s = msg.str();
        throw s.c_str();
    }
    if( align < 0 || align >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "align=" << align << ", modulus=" << modulus << std::endl;
        const std::string s = msg.str();
        throw s.c_str();
    }
    if( modulus <= 0 )
        throw "Modulus must be positive.";
    PopCallStack();
#endif
    return (index + modulus - align) % modulus;
}

#endif /* ELEMENTAL_UTILITIES_HPP */

