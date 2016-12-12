/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_WRITE_BINARY_HPP
#define EL_WRITE_BINARY_HPP

namespace El {
namespace write {

template<typename T>
inline void
Binary( const Matrix<T>& A, string basename="matrix" )
{
    EL_DEBUG_CSE
    
    string filename = basename + "." + FileExtension(BINARY);
    ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int n;
    n = A.Height();
    file.write( (char*)&n, sizeof(Int) );
    n = A.Width();
    file.write( (char*)&n, sizeof(Int) );
    if( A.Height() == A.LDim() )
        file.write( (char*)A.LockedBuffer(), A.Height()*A.Width()*sizeof(T) );
    else
        for( Int j=0; j<A.Width(); ++j )
            file.write( (char*)A.LockedBuffer(0,j), A.Height()*sizeof(T) );
}

} // namespace write
} // namespace El

#endif // ifndef EL_WRITE_BINARY_HPP
