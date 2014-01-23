/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_READ_HPP
#define ELEM_IO_READ_HPP

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
    LogicError("Not yet written");
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

template<typename T>
inline void
Binary( Matrix<T>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file >> height;
    file >> width;
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


template<typename T,Distribution U,Distribution V>
inline void
Binary( DistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file >> height;
    file >> width;
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
    const Int colShift = A.ColShift(); 
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            const Int localIndex = i+j*height;
            const std::streamoff pos = metaBytes + localIndex*sizeof(T);
            file.seekg( pos );
            file.read( (char*)A.Buffer(iLoc,jLoc), sizeof(T) );
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
BinaryFlat
( DistMatrix<T,U,V>& A, Int height, Int width, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
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
    const Int colShift = A.ColShift(); 
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            const Int localIndex = i+j*height;
            const std::streamoff pos = localIndex*sizeof(T);
            file.seekg( pos );
            file.read( (char*)A.Buffer(iLoc,jLoc), sizeof(T) );
        }
    }
}

template<typename T,Distribution V>
inline void
Binary( DistMatrix<T,STAR,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
    std::ifstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    Int height, width;
    file >> height;
    file >> width;
    const Int numBytes = FileSize( file );
    const Int metaBytes = 2*sizeof(Int);
    const Int dataBytes = height*width*sizeof(T);
    const Int numBytesExp = metaBytes + dataBytes;
    if( numBytes != numBytesExp )
        RuntimeError
        ("Expected file to be ",numBytesExp," bytes but found ",numBytes);

    A.Resize( height, width );
    const Int localWidth = A.LocalWidth();
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        const Int localIndex = j*height;
        const std::streamoff pos = metaBytes + localIndex*sizeof(T);
        file.seekg( pos );
        file.read( (char*)A.Buffer(0,jLoc), height*sizeof(T) );
    }
}

template<typename T,Distribution V>
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
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
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
    file >> height;
    file >> width;
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
    file >> height;
    file >> width;
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
    DEBUG_ONLY(CallStackEntry cse("read::Binary"))
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
    default:
        LogicError("Format unsupported for reading");
    }
}

template<typename T,Distribution U,Distribution V>
inline void
Read
( DistMatrix<T,U,V>& A, const std::string filename, FileFormat format=AUTO,
  bool sequential=false )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename ); 

    if( sequential )
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Read( A_CIRC_CIRC.Matrix(), filename, format );
        A = A_CIRC_CIRC;
    }
    else
    {
        switch( format )
        {
        case BINARY:
            read::Binary( A, filename );
            break;
        default:
            LogicError("Unsupported distributed read format"); 
        }
    }
}

// If requesting [* ,* ] or [o ,o ] distributions, no copy is needed

template<typename T>
inline void
Read
( DistMatrix<T,STAR,STAR>& A, const std::string filename, 
  FileFormat format=AUTO )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename );
    if( A.Grid().VCRank() == 0 )
        Read( A.Matrix(), filename, format );
}
template<typename T>
inline void
Read
( DistMatrix<T,CIRC,CIRC>& A, const std::string filename, 
  FileFormat format=AUTO )
{
    DEBUG_ONLY(CallStackEntry cse("Read"))
    if( format == AUTO )
        format = DetectFormat( filename );
    if( A.CrossRank() == A.Root() )
        Read( A.Matrix(), filename, format );
}

} // namespace elem

#endif // ifndef ELEM_IO_READ_HPP
