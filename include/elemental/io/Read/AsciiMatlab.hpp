/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_READ_ASCIIMATLAB_HPP
#define ELEM_READ_ASCIIMATLAB_HPP

namespace elem {
namespace read {

template<typename T>
inline void
AsciiMatlab( Matrix<T>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::AsciiMatlab"))
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);
    LogicError("Not yet written");
}

template<typename T,Dist U,Dist V>
inline void
AsciiMatlab( DistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::AsciiMatlab"))
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);
    LogicError("Not yet written");
}

template<typename T,Dist U,Dist V>
inline void
AsciiMatlab( BlockDistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::AsciiMatlab"))
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);
    LogicError("Not yet written");
}

} // namespace read
} // namespace elem

#endif // ifndef ELEM_READ_ASCIIMATLAB_HPP
