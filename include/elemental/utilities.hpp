/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE.
*/
#ifndef ELEMENTAL_UTILITIES_HPP
#define ELEMENTAL_UTILITIES_HPP 1

namespace elemental {
namespace utilities {

int
GCD
( int a, int b ); 

int
LocalLength
( int n, int shift, int modulus );

int
MaxLocalLength
( int n, int modulus );

int
Shift
( int index, int align, int modulus );

} // utilities
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline int
elemental::utilities::GCD
( int a, int b )
{
#ifndef RELEASE
    if( a < 0 || b < 0 )
        throw std::logic_error( "GCD called with negative argument." );
#endif
    if( b == 0 )
        return a;
    else
        return GCD( b, a-b*(a/b) );
}

inline int
elemental::utilities::LocalLength
( int n, int shift, int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::LocalLength");
    if( n < 0 )
        throw std::logic_error( "n must be non-negative." );
    if( shift < 0 || shift >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid shift: "
            << "shift=" << shift << ", modulus=" << modulus << std::endl;
        throw std::logic_error( msg.str() );
    }
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return ( n > shift ? (n - shift - 1)/modulus + 1 : 0 );
}

inline int
elemental::utilities::MaxLocalLength
( int n, int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::MaxLocalLength");
    if( n < 0 )
        throw std::logic_error( "n must be non-negative." );
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return ( n > 0 ? (n - 1)/modulus + 1 : 0 );
}

// For determining the first global element of process row/column 
// 'index', with distribution alignment 'align' and number of process 
// rows/cols 'modulus'
inline int
elemental::utilities::Shift
( int index, int align, int modulus )
{
#ifndef RELEASE
    PushCallStack("utilities::Shift");
    if( index < 0 || index >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid index: "
            << "index=" << index << ", modulus=" << modulus << std::endl;
        throw std::logic_error( msg.str() );
    }
    if( align < 0 || align >= modulus )
    {
        std::ostringstream msg;
        msg << "Invalid alignment: "
            << "align=" << align << ", modulus=" << modulus << std::endl;
        throw std::logic_error( msg.str() );
    }
    if( modulus <= 0 )
        throw std::logic_error( "Modulus must be positive." );
    PopCallStack();
#endif
    return (index + modulus - align) % modulus;
}

#endif /* ELEMENTAL_UTILITIES_HPP */

