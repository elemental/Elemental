/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef IO_PRINT_HPP
#define IO_PRINT_HPP

namespace elem {

template<typename T>
inline void
Print( const Matrix<T>& A, std::string title="", std::ostream& os=std::cout )
{
#ifndef RELEASE
    CallStackEntry entry("Print");
#endif
    if( title != "" )
        os << title << std::endl;
    
    const int height = A.Height();
    const int width = A.Width();
    for( int i=0; i<height; ++i )
    {
        for( int j=0; j<width; ++j )
            os << A.Get(i,j) << " ";
        os << std::endl;
    }
    os << std::endl;
}

template<typename T,Distribution U,Distribution V>
inline void
Print
( const DistMatrix<T,U,V>& A, std::string title="", std::ostream& os=std::cout )
{
#ifndef RELEASE
    CallStackEntry entry("Print"); 
#endif
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
    if( A.Grid().VCRank() == A_CIRC_CIRC.Root() )
        Print( A_CIRC_CIRC.LockedMatrix(), title, os );
}

// If already in [* ,* ] or [o ,o ] distributions, no copy is needed
template<typename T>
inline void
Print
( const DistMatrix<T,STAR,STAR>& A, std::string title="", 
  std::ostream& os=std::cout )
{
#ifndef RELEASE
    CallStackEntry entry("Print"); 
#endif
    if( A.Grid().VCRank() == 0 )
        Print( A.LockedMatrix(), title, os );
}
template<typename T>
inline void
Print
( const DistMatrix<T,CIRC,CIRC>& A, std::string title="", 
  std::ostream& os=std::cout )
{
#ifndef RELEASE
    CallStackEntry entry("Print"); 
#endif
    if( A.Grid().VCRank() == A.Root() )
        Print( A.LockedMatrix(), title, os );
}

} // namespace elem

#endif // ifndef IO_PRINT_HPP
