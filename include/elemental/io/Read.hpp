/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_READ_HPP
#define ELEM_READ_HPP

namespace elem {
namespace read {

template<typename T>
inline void
Ascii( Matrix<T>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Ascii"))
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Walk through the file once to both count the number of rows and
    // columns and to ensure that the number of columns is consistent
    Int height=0, width=0;
    std::string line;
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

template<typename T,Dist U,Dist V>
inline void
Ascii( DistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Ascii"))
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Walk through the file once to both count the number of rows and
    // columns and to ensure that the number of columns is consistent
    Int height=0, width=0;
    std::string line;
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

template<typename T>
inline void
Binary( Matrix<T>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file.read( (char*)&height, sizeof(Int) );
    file.read( (char*)&width,  sizeof(Int) );
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
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
Binary( DistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file.read( (char*)&height, sizeof(Int) );
    file.read( (char*)&width,  sizeof(Int) );
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            const Int localIndex = i+j*height;
            const std::streamoff pos = metaBytes + localIndex*sizeof(T);
            file.seekg( pos );
            file.read( (char*)A.Buffer(iLoc,jLoc), sizeof(T) );
        }
    }
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

template<typename T,Dist V>
inline void
Binary( DistMatrix<T,STAR,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file.read( (char*)&height, sizeof(Int) );
    file.read( (char*)&width,  sizeof(Int) );
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int localIndex = j*height;
        const std::streamoff pos = metaBytes + localIndex*sizeof(T);
        file.seekg( pos );
        file.read( (char*)A.Buffer(0,jLoc), height*sizeof(T) );
    }
}

template<typename T,Dist V>
inline void
BinaryFlat
( DistMatrix<T,STAR,V>& A, Int height, Int width, const std::string filename )
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

template<typename T>
inline void
Binary( DistMatrix<T,STAR,STAR>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file.read( (char*)&height, sizeof(Int) );
    file.read( (char*)&width,  sizeof(Int) );
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
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

template<typename T>
inline void
BinaryFlat
( DistMatrix<T,STAR,STAR>& A, Int height, Int width, 
  const std::string filename )
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

template<typename T>
inline void
Binary( DistMatrix<T,CIRC,CIRC>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file.read( (char*)&height, sizeof(Int) );
    file.read( (char*)&width,  sizeof(Int) );
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    if( A.CrossRank() == A.Root() )
    {
        if( A.Height() == A.LDim() )
            file.read( (char*)A.Buffer(), height*width*sizeof(T) );
        else
            for( Int j=0; j<width; ++j )
                file.read( (char*)A.Buffer(0,j), height*sizeof(T) );
    }
}

template<typename T>
inline void
BinaryFlat
( DistMatrix<T,CIRC,CIRC>& A, Int height, Int width, 
  const std::string filename )
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
    if( A.CrossRank() == A.Root() )
    {
        if( A.Height() == A.LDim() )
            file.read( (char*)A.Buffer(), height*width*sizeof(T) );
        else
            for( Int j=0; j<width; ++j )
                file.read( (char*)A.Buffer(0,j), height*sizeof(T) );
    }
}

} // namespace read

template<typename T>
inline void
Read( Matrix<T>& A, const std::string filename, FileFormat format=AUTO )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename );

    switch( format )
    {
    case ASCII:
        read::Ascii( A, filename );
        break;
    case ASCII_MATLAB:
        read::AsciiMatlab( A, filename );
        break;
    case BINARY:
        read::Binary( A, filename );
        break;
    case BINARY_FLAT:
        read::BinaryFlat( A, A.Height(), A.Width(), filename );
        break;
    default:
        LogicError("Format unsupported for reading");
    }
}

template<typename T,Dist U,Dist V>
inline void
Read
( DistMatrix<T,U,V>& A, const std::string filename, FileFormat format=AUTO,
  bool sequential=false )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename ); 

    if( U == A.UGath && V == A.VGath )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Read( A.Matrix(), filename, format );
    }
    else if( sequential )
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
        if( format == BINARY_FLAT )
            A_CIRC_CIRC.Resize( A.Height(), A.Width() );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Read( A_CIRC_CIRC.Matrix(), filename, format );
        A = A_CIRC_CIRC;
    }
    else
    {
        switch( format )
        {
        case ASCII:
            read::Ascii( A, filename );
            break;
        case ASCII_MATLAB:
            read::AsciiMatlab( A, filename );
            break;
        case BINARY:
            read::Binary( A, filename );
            break;
        case BINARY_FLAT:
            read::BinaryFlat( A, A.Height(), A.Width(), filename );
            break;
        default:
            LogicError("Unsupported distributed read format"); 
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_READ_HPP
