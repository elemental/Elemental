/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_READ_ASCIIMATLAB_HPP
#define EL_READ_ASCIIMATLAB_HPP

namespace El {
namespace read {

template<typename T>
inline void
AsciiMatlab( Matrix<T>& A, const string filename )
{
    EL_DEBUG_CSE
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);
    LogicError("Not yet written");
}

template<typename T>
inline void
AsciiMatlab( AbstractDistMatrix<T>& A, const string filename )
{
    EL_DEBUG_CSE
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);
    LogicError("Not yet written");
}

} // namespace read
} // namespace El

#endif // ifndef EL_READ_ASCIIMATLAB_HPP
