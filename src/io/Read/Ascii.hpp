/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_READ_ASCII_HPP
#define EL_READ_ASCII_HPP

namespace El {
namespace read {

template<typename T>
inline void
Ascii( Matrix<T>& A, const string filename )
{
    EL_DEBUG_CSE
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Walk through the file once to both count the number of rows and
    // columns and to ensure that the number of columns is consistent
    Int height=0, width=0;
    string line;
    while( std::getline( file, line ) )
    {
        std::stringstream lineStream( line );
        Int numCols=0;
        T value;
        while( lineStream >> value ) ++numCols;
        if( numCols != 0 )
        {
            if( numCols != width && width != 0 )
                LogicError("Inconsistent number of columns");
            else
                width = numCols;
            ++height;
        }
    }
    file.clear();
    file.seekg(0,file.beg);

    // Resize the matrix and then read it
    A.Resize( height, width );
    Int i=0;
    while( std::getline( file, line ) )
    {
        std::stringstream lineStream( line );
        Int j=0;
        T value;
        while( lineStream >> value )
        {
            A.Set( i, j, value );
            ++j;
        }
        ++i;
    }
}

template<typename T>
inline void
Ascii( AbstractDistMatrix<T>& A, const string filename )
{
    EL_DEBUG_CSE
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Walk through the file once to both count the number of rows and
    // columns and to ensure that the number of columns is consistent
    Int height=0, width=0;
    string line;
    while( std::getline( file, line ) )
    {
        std::stringstream lineStream( line );
        Int numCols=0;
        T value;
        while( lineStream >> value ) ++numCols;
        if( numCols != 0 )
        {
            if( numCols != width && width != 0 )
                LogicError("Inconsistent number of columns");
            else
                width = numCols;
            ++height;
        }
    }
    file.clear();
    file.seekg(0,file.beg);

    // Resize the matrix and then read in our local portion
    A.Resize( height, width );
    Int i=0;
    while( std::getline( file, line ) )
    {
        std::stringstream lineStream( line );
        Int j=0;
        T value;
        while( lineStream >> value )
        {
            A.Set( i, j, value );
            ++j;
        }
        ++i;
    }
}

} // namespace read
} // namespace El

#endif // ifndef EL_READ_ASCII_HPP
