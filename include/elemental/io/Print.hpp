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
    // TODO: Optimize to only gather to the root
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    if( A.Grid().VCRank() == 0 )
        Print( A_STAR_STAR.LockedMatrix(), title, os );
}

} // namespace elem

#endif // ifndef IO_PRINT_HPP
