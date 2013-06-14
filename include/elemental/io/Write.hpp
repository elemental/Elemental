/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef IO_WRITE_HPP
#define IO_WRITE_HPP

#include "elemental/io/Print.hpp"

namespace elem {

template<typename T>
inline void
Write( const Matrix<T>& A, std::string title="", std::string filename="Matrix" )
{
#ifndef RELEASE
    CallStackEntry entry("Write");
#endif
    std::ofstream file( filename.c_str() );
    file.setf( std::ios::scientific );
    Print( A, title, file );
    file.close();
}

template<typename T,Distribution U,Distribution V>
inline void
Write
( const DistMatrix<T,U,V>& A, std::string title="", 
  std::string filename="DistMatrix" )
{
#ifndef RELEASE
    CallStackEntry entry("Write"); 
#endif
    // TODO: Optimize to only gather to the root
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    if( A.Grid().VCRank() == 0 )
        Write( A_STAR_STAR.LockedMatrix(), title, filename );
}

} // namespace elem

#endif // ifndef IO_WRITE_HPP
