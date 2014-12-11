/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENVIRONMENT_IMPL_HPP
#define EL_ENVIRONMENT_IMPL_HPP

namespace El {

inline void Args::HandleVersion( std::ostream& os ) const
{
    std::string version = "--version";
    char** arg = std::find( argv_, argv_+argc_, version );
    const bool foundVersion = ( arg != argv_+argc_ );
    if( foundVersion )
    {
        if( mpi::WorldRank() == 0 )
            PrintVersion();
        throw ArgException();
    }
}

inline void Args::HandleBuild( std::ostream& os ) const
{
    std::string build = "--build";
    char** arg = std::find( argv_, argv_+argc_, build );
    const bool foundBuild = ( arg != argv_+argc_ );
    if( foundBuild )
    {
        if( mpi::WorldRank() == 0 )
        {
            PrintVersion();
            PrintConfig();
            PrintCCompilerInfo();
            PrintCxxCompilerInfo();
        }
        throw ArgException();
    }
}

template<typename T>
inline T
Input( std::string name, std::string desc )
{ return GetArgs().Input<T>( name, desc ); }

template<typename T>
inline T
Input( std::string name, std::string desc, T defaultVal )
{ return GetArgs().Input( name, desc, defaultVal ); }

inline void
ProcessInput()
{ GetArgs().Process(); }

inline void
PrintInputReport()
{ GetArgs().PrintReport(); }

inline void ReportException( const std::exception& e, std::ostream& os )
{
    try {
        const ArgException& argExcept = dynamic_cast<const ArgException&>(e);
        if( std::string(argExcept.what()) != "" ) 
            os << argExcept.what() << std::endl;
        DEBUG_ONLY(DumpCallStack(os))
    } 
    catch( std::exception& castExcept ) 
    { 
        if( std::string(e.what()) != "" )
        {
            os << "Process " << mpi::WorldRank() << " caught error message:\n"
               << e.what() << std::endl;
        }
        DEBUG_ONLY(DumpCallStack(os))
        mpi::Abort( mpi::COMM_WORLD, 1 );
    }
}

inline void ComplainIfDebug()
{
    DEBUG_ONLY(
        if( mpi::WorldRank() == 0 )
        {
            std::cout << "==========================================\n"
                      << " In debug mode! Performance will be poor! \n"
                      << "==========================================" 
                      << std::endl;
        }
    )
}

template<typename T>
inline void 
MemCopy( T* dest, const T* source, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memcpy( dest, source, numEntries*sizeof(T) );
}

template<typename T>
inline void
MemSwap( T* a, T* b, T* temp, std::size_t numEntries )
{
    // temp := a
    std::memcpy( temp, a, numEntries*sizeof(T) );
    // a := b
    std::memcpy( a, b, numEntries*sizeof(T) );
    // b := temp
    std::memcpy( b, temp, numEntries*sizeof(T) );
}

template<typename T>
inline void
StridedMemCopy
(       T* dest,   Int destStride, 
  const T* source, Int sourceStride, Int numEntries )
{
    // For now, use the BLAS wrappers/generalization
    blas::Copy( numEntries, source, sourceStride, dest, destStride );
}

template<typename T>
inline void
InterleaveMatrix
( Int localHeight, Int localWidth,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB )
{
    // TODO: Add OpenMP parallelization and/or optimize
    if( colStrideA == 1 && colStrideB == 1 )
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            MemCopy( &B[jLoc*rowStrideB], &A[jLoc*rowStrideA], localHeight );
    }
    else
    {
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            StridedMemCopy
            ( &B[jLoc*rowStrideB], colStrideB,
              &A[jLoc*rowStrideA], colStrideA, localHeight );
    }
}

template<typename T>
inline void
InterleaveMatrixUpdate
( T alpha, Int localHeight, Int localWidth,
  const T* A, Int colStrideA, Int rowStrideA,
        T* B, Int colStrideB, Int rowStrideB )
{
    // TODO: Add OpenMP parallelization and/or optimize
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        blas::Axpy
        ( localHeight, alpha, 
          &A[rowStrideA*jLoc], colStrideA,
          &B[rowStrideB*jLoc], rowStrideB );
}

template<typename T>
inline void 
MemZero( T* buffer, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memset( buffer, 0, numEntries*sizeof(T) );
}

template<typename T>
inline void SwapClear( T& x ) { T().swap( x ); }

template<typename T>
inline T 
Scan( const std::vector<T>& counts, std::vector<T>& offsets )
{
    offsets.resize( counts.size() );
    T total = 0;
    for( Int i=0; i<counts.size(); ++i )
    {
        offsets[i] = total;
        total += counts[i];
    }
    return total;
}

template<typename T>
inline void
EnsureConsistent( T alpha, mpi::Comm comm, std::string name )
{
    std::string tag = ( name=="" ? "" : name+" " );
    const Int commSize = mpi::Size( comm );
    const Int commRank = mpi::Rank( comm );
    std::vector<T> a(commSize);
    mpi::Gather( &alpha, 1, a.data(), 1, 0, comm );
    if( commRank == 0 ) 
    {
        for( Int j=0; j<commSize; ++j )
            if( a[j] != alpha )
                std::cout << "Process " << j << "'s " << tag << "value, " 
                          << a[j] << ", mismatched the root's, " << alpha 
                          << std::endl;
    }
}

template<typename F>
inline void UpdateScaledSquare( F alpha, Base<F>& scale, Base<F>& scaledSquare )
{
    typedef Base<F> Real;
    Real alphaAbs = Abs(alpha);
    if( alphaAbs != 0 )
    {
        if( alphaAbs <= scale )
        {
            const Real relScale = alphaAbs/scale;
            scaledSquare += relScale*relScale;
        }
        else
        {
            const Real relScale = scale/alphaAbs;
            scaledSquare = scaledSquare*relScale*relScale + Real(1);
            scale = alphaAbs;
        }
    }
}

} // namespace El

#endif // ifndef EL_ENVIRONMENT_IMPL_HPP
