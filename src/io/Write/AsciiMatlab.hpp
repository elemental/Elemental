/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_WRITE_ASCIIMATLAB_HPP
#define EL_WRITE_ASCIIMATLAB_HPP

namespace El {
namespace write {

template<typename T>
inline void
AsciiMatlab
( const Matrix<T>& A, string basename="matrix", string title="matrix" )
{
    DEBUG_CSE
    // Empty titles are not legal
    if( title == "" )
        title = "matrix";

    string filename = basename + "." + FileExtension(ASCII_MATLAB);
    ofstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    file.setf( std::ios::scientific );
    file << title << " = [\n";
    Print( A, "", file );
    file << "];\n";
}

} // namespace write
} // namespace El

#endif // ifndef EL_WRITE_ASCIIMATLAB_HPP
