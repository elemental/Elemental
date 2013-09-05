/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_ENVIRONMENT_IMPL_HPP
#define ELEM_CORE_ENVIRONMENT_IMPL_HPP

namespace elem {

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
    if( std::string(e.what()) != "" )
    {
        os << "Process " << mpi::WorldRank() << " caught error message:\n"
           << e.what() << std::endl;
    }
#ifndef RELEASE
    DumpCallStack( os );
#endif
}

inline void ComplainIfDebug()
{
#ifndef RELEASE
    if( mpi::WorldRank() == 0 )
    {
        std::cout << "==========================================\n"
                  << " In debug mode! Performance will be poor! \n"
                  << "==========================================" << std::endl;
    }
#endif
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
(       T* dest,   std::size_t destStride, 
  const T* source, std::size_t sourceStride, std::size_t numEntries )
{
    // For now, use the BLAS wrappers/generalization
    blas::Copy( numEntries, source, sourceStride, dest, destStride );
}

template<typename T>
inline void 
MemZero( T* buffer, std::size_t numEntries )
{
    // This can be optimized/generalized later
    std::memset( buffer, 0, numEntries*sizeof(T) );
}

} // namespace elem

#endif // ifndef ELEM_CORE_ENVIRONMENT_IMPL_HPP
