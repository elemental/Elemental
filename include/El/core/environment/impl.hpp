/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ENVIRONMENT_IMPL_HPP
#define EL_ENVIRONMENT_IMPL_HPP

namespace El {

template<typename T>
inline T
Input( string name, string desc )
{ return GetArgs().Input<T>( name, desc ); }

template<typename T>
inline T
Input( string name, string desc, T defaultVal )
{ return GetArgs().Input( name, desc, defaultVal ); }

inline void
ProcessInput()
{ GetArgs().Process(); }

inline void
PrintInputReport()
{ GetArgs().PrintReport(); }

template<typename T>
inline void 
MemCopy( T* dest, const T* source, size_t numEntries )
{
    // This can be optimized/generalized later
    std::memcpy( dest, source, numEntries*sizeof(T) );
}
#ifdef EL_HAVE_MPC
inline void
MemCopy( BigInt* dest, const BigInt* source, size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        dest[k] = source[k];
}

inline void
MemCopy( BigFloat* dest, const BigFloat* source, size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        dest[k] = source[k];
}
#endif

template<typename T>
inline void
MemSwap( T* a, T* b, T* temp, size_t numEntries )
{
    // temp := a
    MemCopy( temp, a, numEntries );
    // a := b
    MemCopy( a, b, numEntries );
    // b := temp
    MemCopy( b, temp, numEntries );
}
#ifdef EL_HAVE_MPC
inline void
MemSwap( BigInt* a, BigInt* b, BigInt* temp, size_t numEntries )
{
    // NOTE: This is the same as above for now

    // temp := a
    MemCopy( temp, a, numEntries );
    // a := b
    MemCopy( a, b, numEntries );
    // b := temp
    MemCopy( b, temp, numEntries );
}

inline void
MemSwap( BigFloat* a, BigFloat* b, BigFloat* temp, size_t numEntries )
{
    // NOTE: This is the same as above for now

    // temp := a
    MemCopy( temp, a, numEntries );
    // a := b
    MemCopy( a, b, numEntries );
    // b := temp
    MemCopy( b, temp, numEntries );
}
#endif

template<typename T>
inline void
StridedMemCopy
(       T* dest,   Int destStride, 
  const T* source, Int sourceStride, Int numEntries )
{
    // For now, use the BLAS wrappers/generalization
    blas::Copy( numEntries, source, sourceStride, dest, destStride );
}
#ifdef EL_HAVE_MPC
inline void
StridedMemCopy
(       BigInt* dest,   Int destStride,
  const BigInt* source, Int sourceStride, Int numEntries )
{
    for( Int k=0; k<numEntries; ++k )
        dest[destStride*k] = source[sourceStride*k];
}

inline void
StridedMemCopy
(       BigFloat* dest,   Int destStride,
  const BigFloat* source, Int sourceStride, Int numEntries )
{
    for( Int k=0; k<numEntries; ++k )
        dest[destStride*k] = source[sourceStride*k];
}
#endif

template<typename T>
inline void 
MemZero( T* buffer, size_t numEntries )
{
    // This can be optimized/generalized later
    std::memset( buffer, 0, numEntries*sizeof(T) );
}
#ifdef EL_HAVE_MPC
inline void MemZero( BigInt* buffer, size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        buffer[k].Zero();
}

inline void MemZero( BigFloat* buffer, size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        buffer[k].Zero();
}
#endif

template<typename T>
inline void SwapClear( T& x ) { T().swap( x ); }

template<typename T>
inline T 
Scan( const vector<T>& counts, vector<T>& offsets )
{
    offsets.resize( counts.size() );
    T total = 0;
    for( size_t i=0; i<counts.size(); ++i )
    {
        offsets[i] = total;
        total += counts[i];
    }
    return total;
}

template<typename T>
inline void
EnsureConsistent( T alpha, mpi::Comm comm, string name )
{
    string tag = ( name=="" ? "" : name+" " );
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    vector<T> a(commSize);
    mpi::Gather( &alpha, 1, a.data(), 1, 0, comm );
    if( commRank == 0 ) 
    {
        for( Int j=0; j<commSize; ++j )
            if( a[j] != alpha )
                cout << "Process " << j << "'s " << tag << "value, " 
                     << a[j] << ", mismatched the root's, " << alpha 
                     << endl;
    }
}

} // namespace El

#endif // ifndef EL_ENVIRONMENT_IMPL_HPP
