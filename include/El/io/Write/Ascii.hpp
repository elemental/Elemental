/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_WRITE_ASCII_HPP
#define EL_WRITE_ASCII_HPP

#include "../Print.hpp"

namespace El {
namespace write {

template<typename T>
inline void
Ascii( const Matrix<T>& A, std::string basename="matrix", std::string title="" )
{
    DEBUG_ONLY(CallStackEntry cse("write::Ascii"))
    std::string filename = basename + "." + FileExtension(ASCII);
    std::ofstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    file.setf( std::ios::scientific );
    Print( A, title, file );
}

} // namespace write
} // namespace El

#endif // ifndef EL_WRITE_ASCII_HPP
