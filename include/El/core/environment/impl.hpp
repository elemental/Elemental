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
T Input( string name, string desc )
{ return GetArgs().Input<T>( name, desc ); }

template<typename T>
T Input( string name, string desc, T defaultVal )
{ return GetArgs().Input( name, desc, defaultVal ); }

inline void
ProcessInput()
{ GetArgs().Process(); }

inline void
PrintInputReport()
{ GetArgs().PrintReport(); }

template<typename T,
         typename/*=EnableIf<IsScalar<T>>*/>
const T& Max( const T& m, const T& n ) EL_NO_EXCEPT
{ return std::max(m,n); }

inline const Int& Max( const Int& m, const Int& n ) EL_NO_EXCEPT
{ return std::max(m,n); }

template<typename T,
         typename/*=EnableIf<IsScalar<T>>*/>
const T& Min( const T& m, const T& n ) EL_NO_EXCEPT
{ return std::min(m,n); }

inline const Int& Min( const Int& m, const Int& n ) EL_NO_EXCEPT
{ return std::min(m,n); }

template<typename T,
         typename/*=EnableIf<IsPacked<T>>*/>
void MemCopy
(       T* dest,
  const T* source,
        size_t numEntries )
{
    // This can be optimized/generalized later
    std::memcpy( dest, source, numEntries*sizeof(T) );
}
template<typename T,
         typename/*=DisableIf<IsPacked<T>>*/,
         typename/*=void*/>
void MemCopy
(       T* dest,
  const T* source,
        size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        dest[k] = source[k];
}

template<typename T,
         typename/*=EnableIf<IsPacked<T>>*/>
void MemSwap( T* a, T* b, T* temp, size_t numEntries )
{
    // temp := a
    MemCopy( temp, a, numEntries );
    // a := b
    MemCopy( a, b, numEntries );
    // b := temp
    MemCopy( b, temp, numEntries );
}
template<typename T,
         typename/*=DisableIf<IsPacked<T>>*/,
         typename/*=void*/>
void MemSwap
( T* a,
  T* b,
  T* temp,
  size_t numEntries )
{
    // TODO: Optimize
    // temp := a
    MemCopy( temp, a, numEntries );
    // a := b
    MemCopy( a, b, numEntries );
    // b := temp
    MemCopy( b, temp, numEntries );
}

template<typename T,
         typename/*=EnableIf<IsPacked<T>>*/>
void StridedMemCopy
(       T* dest,   Int destStride,
  const T* source, Int sourceStride, Int numEntries )
{
    // For now, use the BLAS wrappers/generalization
    blas::Copy( numEntries, source, sourceStride, dest, destStride );
}
template<typename T,
         typename/*=DisableIf<IsPacked<T>>*/,
         typename/*=void*/>
void StridedMemCopy
(       T* dest,   Int destStride,
  const T* source, Int sourceStride, Int numEntries )
{
    for( Int k=0; k<numEntries; ++k )
        dest[destStride*k] = source[sourceStride*k];
}

template<typename S,typename T>
void CopySTL( const S& a, T& b )
{
    b.resize( a.size() );
    std::copy( a.begin(), a.end(), b.begin() );
}

template<typename T,
         typename/*=EnableIf<IsPacked<T>>*/>
void MemZero( T* buffer, size_t numEntries )
{
    // This can be optimized/generalized later
    std::memset( buffer, 0, numEntries*sizeof(T) );
}
template<typename T,
         typename/*=DisableIf<IsPacked<T>>*/,
         typename/*=void*/>
void MemZero( T* buffer, size_t numEntries )
{
    for( size_t k=0; k<numEntries; ++k )
        buffer[k].Zero();
}

template<typename T>
void SwapClear( T& x ) { T().swap( x ); }

template<typename T,
         typename/*=EnableIf<IsPacked<T>>*/>
void FastResize( vector<T>& v, Int numEntries )
{
#ifdef EL_ZERO_INIT
    v.resize( numEntries );
#elif defined(EL_HAVE_VALGRIND)
    if( EL_RUNNING_ON_VALGRIND )
        v.resize( numEntries );
    else
        v.reserve( numEntries );
#else
    v.reserve( numEntries );
#endif
}
template<typename T,
         typename/*=DisableIf<IsPacked<T>>*/,
         typename/*=void*/>
void FastResize( vector<T>& v, Int numEntries )
{ v.resize( numEntries ); }

template<typename T,typename... ArgPack>
void BuildStream( ostringstream& os, const T& item, const ArgPack& ... args )
{
    os << item;
    BuildStream( os, args... );
}

template<typename... ArgPack>
string BuildString( const ArgPack& ... args )
{
    ostringstream os;
    BuildStream( os, args... );
    return os.str();
}

template<typename... ArgPack>
void UnrecoverableError( const ArgPack& ... args )
{
    ostringstream os;
    BuildStream( os, args... );
    os << endl;
    UnrecoverableException( os.str().c_str() );
}

template<typename... ArgPack>
void LogicError( const ArgPack& ... args )
{
    ostringstream os;
    BuildStream( os, args... );
    os << endl;
    throw std::logic_error( os.str().c_str() );
}

template<typename... ArgPack>
void RuntimeError( const ArgPack& ... args )
{
    ostringstream os;
    BuildStream( os, args... );
    os << endl;
    throw std::runtime_error( os.str().c_str() );
}

template<class MatType>
string DimsString( const MatType& A, string label )
{
    ostringstream os;
    os << label << " ~ " << A.Height() << " x " << A.Width();
    return os.str();
}

template<typename... ArgPack>
void Log( const ArgPack& ... args )
{
    std::ostringstream str;
    BuildStream( str, args... );
    LogOS() << str.str() << std::endl;
}

template<typename... ArgPack>
void Output( const ArgPack& ... args )
{
    ostringstream os;
    os << Indent();
    BuildStream( os, args... );
    os << endl;
    cout << os.str();
}

template<typename... ArgPack>
void OutputFromRoot( mpi::Comm comm, const ArgPack& ... args )
{
    if( mpi::Rank(comm) == 0 )
    {
        Output( args... );
    }
}

template<typename T>
T Scan( const vector<T>& counts, vector<T>& offsets )
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
void EnsureConsistent( T alpha, mpi::Comm comm, string name )
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
