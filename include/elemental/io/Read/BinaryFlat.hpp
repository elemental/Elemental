/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_READ_BINARYFLAT_HPP
#define ELEM_READ_BINARYFLAT_HPP

namespace elem {
namespace read {

template<typename T>
inline void
BinaryFlat( Matrix<T>& A, Int height, Int width, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::BinaryFlat"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    const Int numBytes = FileSize( file );
    const Int numBytesExp = height*width*sizeof(T);
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    if( A.Height() == A.LDim() )
        file.read( (char*)A.Buffer(), height*width*sizeof(T) );
    else
        for( Int j=0; j<width; ++j )
            file.read( (char*)A.Buffer(0,j), height*sizeof(T) );
}

template<typename T,Dist U,Dist V>
inline void
BinaryFlat
( DistMatrix<T,U,V>& A, Int height, Int width, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::BinaryFlat"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    const Int numBytes = FileSize( file );
    const Int numBytesExp = height*width*sizeof(T);
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() )
        {
            if( A.Height() == A.LDim() )
                file.read( (char*)A.Buffer(), height*width*sizeof(T) );
            else
                for( Int j=0; j<width; ++j )
                    file.read( (char*)A.Buffer(0,j), height*sizeof(T) );
        }
    }
    else if( U == A.UGath )
    {
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int localIndex = j*height;
            const std::streamoff pos = localIndex*sizeof(T);
            file.seekg( pos );
            file.read( (char*)A.Buffer(0,jLoc), height*sizeof(T) );
        }
    }
    else
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int localIndex = i+j*height;
                const std::streamoff pos = localIndex*sizeof(T);
                file.seekg( pos );
                file.read( (char*)A.Buffer(iLoc,jLoc), sizeof(T) );
            }
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
BinaryFlat
( BlockDistMatrix<T,U,V>& A, Int height, Int width, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::BinaryFlat"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    const Int numBytes = FileSize( file );
    const Int numBytesExp = height*width*sizeof(T);
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() )
        {
            if( A.Height() == A.LDim() )
                file.read( (char*)A.Buffer(), height*width*sizeof(T) );
            else
                for( Int j=0; j<width; ++j )
                    file.read( (char*)A.Buffer(0,j), height*sizeof(T) );
        }
    }
    else if( U == A.UGath )
    {
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            const Int localIndex = j*height;
            const std::streamoff pos = localIndex*sizeof(T);
            file.seekg( pos );
            file.read( (char*)A.Buffer(0,jLoc), height*sizeof(T) );
        }
    }
    else
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int localIndex = i+j*height;
                const std::streamoff pos = localIndex*sizeof(T);
                file.seekg( pos );
                file.read( (char*)A.Buffer(iLoc,jLoc), sizeof(T) );
            }
        }
    }
}

} // namespace read
} // namespace elem

#endif // ifndef ELEM_READ_BINARYFLAT_HPP
