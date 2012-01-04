/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"

namespace {
// Core routines
bool initializedElemental = false;
bool elementalInitializedMpi;
std::stack<int> blocksizeStack;
elemental::Grid* defaultGrid = 0;

// Debugging
#ifndef RELEASE
std::stack<std::string> callStack;
#endif

// Tuning paramters for basic routines
int localHemvFloatBlocksize = 64;
int localHemvDoubleBlocksize = 64;
int localHemvComplexFloatBlocksize = 64;
int localHemvComplexDoubleBlocksize = 64;

int localSymvFloatBlocksize = 64;
int localSymvDoubleBlocksize = 64;
int localSymvComplexFloatBlocksize = 64;
int localSymvComplexDoubleBlocksize = 64;

int localTrr2kFloatBlocksize = 64;
int localTrr2kDoubleBlocksize = 64;
int localTrr2kComplexFloatBlocksize = 64;
int localTrr2kComplexDoubleBlocksize = 64;

int localTrrkFloatBlocksize = 64;
int localTrrkDoubleBlocksize = 64;
int localTrrkComplexFloatBlocksize = 64;
int localTrrkComplexDoubleBlocksize = 64;

// Tuning parameters for advanced routines
using namespace elemental;
HermitianTridiagApproach tridiagApproach = HERMITIAN_TRIDIAG_DEFAULT;
GridOrder gridOrder = ROW_MAJOR;
}

namespace elemental {

bool Initialized()
{ return ::initializedElemental; }

void Initialize( int& argc, char**& argv )
{
    // If Elemental is currently initialized, then this is a no-op
    if( ::initializedElemental )
        return;

    if( !mpi::Initialized() )
    {
        if( mpi::Finalized() )
        {
            throw std::logic_error
            ("Cannot initialize elemental after finalizing MPI");
        }
#ifdef _OPENMP
        const int provided = 
            mpi::InitializeThread
            ( argc, argv, mpi::THREAD_MULTIPLE );
        if( provided != mpi::THREAD_MULTIPLE )
        {
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
        }
#else
        mpi::Initialize( argc, argv );
#endif
        ::elementalInitializedMpi = true;
    }
    else
    {
        ::elementalInitializedMpi = false;
    }

    // Queue a default algorithmic blocksize
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::blocksizeStack.push( 128 );

    // Build the default grid
    defaultGrid = new Grid( mpi::COMM_WORLD );

    // Build the pivot operations needed by the distributed LU
    internal::CreatePivotOp<float>();
    internal::CreatePivotOp<double>();
    internal::CreatePivotOp<std::complex<float> >();
    internal::CreatePivotOp<std::complex<double> >();

    // Seed the parallel random number generator, PLCG
    plcg::UInt64 seed;
    seed.d[0] = time(0);
    seed.d[1] = time(0);
    const unsigned rank = mpi::CommRank( mpi::COMM_WORLD );
    const unsigned size = mpi::CommSize( mpi::COMM_WORLD );
    mpi::Broadcast
    ( (byte*)seed.d, 2*sizeof(unsigned), 0, mpi::COMM_WORLD );
    plcg::SeedParallelLcg( rank, size, seed );

    ::initializedElemental = true;
}

void Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    // If Elemental is not currently initialized, then this is a no-op
    if( !::initializedElemental )
        return;

    if( mpi::Finalized() )
    {
        std::cerr << "Warning: MPI was finalized before Elemental." 
                  << std::endl;
    }
    else if( ::elementalInitializedMpi )
    {
        // Destroy the pivot ops needed by the distributed LU
        internal::DestroyPivotOp<float>();
        internal::DestroyPivotOp<double>();
        internal::DestroyPivotOp<std::complex<float> >();
        internal::DestroyPivotOp<std::complex<double> >();

        // Delete the default grid
        delete ::defaultGrid;
        ::defaultGrid = 0;

        mpi::Finalize();
    }

    delete ::defaultGrid;
    ::defaultGrid = 0;
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::initializedElemental = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

int Blocksize()
{ return ::blocksizeStack.top(); }

void SetBlocksize( int blocksize )
{ ::blocksizeStack.top() = blocksize; }

void PushBlocksizeStack( int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void PopBlocksizeStack()
{ ::blocksizeStack.pop(); }

const Grid& DefaultGrid()
{
#ifndef RELEASE
    PushCallStack("DefaultGrid");
    if( ::defaultGrid == 0 )
        throw std::logic_error
        ("Attempted to return a non-existant default grid. Please ensure that "
         "Elemental is initialized before creating a DistMatrix.");
    PopCallStack();
#endif
    return *::defaultGrid;
}

// If we are not in RELEASE mode, then implement wrappers for a CallStack
#ifndef RELEASE
void PushCallStack( std::string s )
{ ::callStack.push(s); }

void PopCallStack()
{ ::callStack.pop(); }

void DumpCallStack()
{
    std::ostringstream msg;
    while( ! ::callStack.empty() )
    {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n";
        ::callStack.pop();
    }
    std::cerr << msg.str() << std::endl;
}
#endif // RELEASE

template<>
void SetLocalHemvBlocksize<float>( int blocksize )
{ ::localHemvFloatBlocksize = blocksize; }

template<>
void SetLocalHemvBlocksize<double>( int blocksize )
{ ::localHemvDoubleBlocksize = blocksize; }

template<>
void SetLocalHemvBlocksize<std::complex<float> >( int blocksize )
{ ::localHemvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalHemvBlocksize<std::complex<double> >( int blocksize )
{ ::localHemvComplexDoubleBlocksize = blocksize; }

template<>
int LocalHemvBlocksize<float>()
{ return ::localHemvFloatBlocksize; }

template<>
int LocalHemvBlocksize<double>()
{ return ::localHemvDoubleBlocksize; }

template<>
int LocalHemvBlocksize<std::complex<float> >()
{ return ::localHemvComplexFloatBlocksize; }

template<>
int LocalHemvBlocksize<std::complex<double> >()
{ return ::localHemvComplexDoubleBlocksize; }

template<>
void SetLocalSymvBlocksize<float>( int blocksize )
{ ::localSymvFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<double>( int blocksize )
{ ::localSymvDoubleBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<std::complex<float> >( int blocksize )
{ ::localSymvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<std::complex<double> >( int blocksize )
{ ::localSymvComplexDoubleBlocksize = blocksize; }

template<>
int LocalSymvBlocksize<float>()
{ return ::localSymvFloatBlocksize; }

template<>
int LocalSymvBlocksize<double>()
{ return ::localSymvDoubleBlocksize; }

template<>
int LocalSymvBlocksize<std::complex<float> >()
{ return ::localSymvComplexFloatBlocksize; }

template<>
int LocalSymvBlocksize<std::complex<double> >()
{ return ::localSymvComplexDoubleBlocksize; }

template<>
void SetLocalTrr2kBlocksize<float>( int blocksize )
{ ::localTrr2kFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<double>( int blocksize )
{ ::localTrr2kDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<std::complex<float> >( int blocksize )
{ ::localTrr2kComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<std::complex<double> >( int blocksize )
{ ::localTrr2kComplexDoubleBlocksize = blocksize; }

template<>
int LocalTrr2kBlocksize<float>()
{ return ::localTrr2kFloatBlocksize; }

template<>
int LocalTrr2kBlocksize<double>()
{ return ::localTrr2kDoubleBlocksize; }

template<>
int LocalTrr2kBlocksize<std::complex<float> >()
{ return ::localTrr2kComplexFloatBlocksize; }

template<>
int LocalTrr2kBlocksize<std::complex<double> >()
{ return ::localTrr2kComplexDoubleBlocksize; }

template<>
void SetLocalTrrkBlocksize<float>( int blocksize )
{ ::localTrrkFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<double>( int blocksize )
{ ::localTrrkDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<std::complex<float> >( int blocksize )
{ ::localTrrkComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<std::complex<double> >( int blocksize )
{ ::localTrrkComplexDoubleBlocksize = blocksize; }

template<>
int LocalTrrkBlocksize<float>()
{ return ::localTrrkFloatBlocksize; }

template<>
int LocalTrrkBlocksize<double>()
{ return ::localTrrkDoubleBlocksize; }

template<>
int LocalTrrkBlocksize<std::complex<float> >()
{ return ::localTrrkComplexFloatBlocksize; }

template<>
int LocalTrrkBlocksize<std::complex<double> >()
{ return ::localTrrkComplexDoubleBlocksize; }

void SetHermitianTridiagApproach( HermitianTridiagApproach approach )
{ ::tridiagApproach = approach; }

HermitianTridiagApproach GetHermitianTridiagApproach()
{ return ::tridiagApproach; }

void SetHermitianTridiagGridOrder( GridOrder order )
{ ::gridOrder = order; }

GridOrder GetHermitianTridiagGridOrder()
{ return ::gridOrder; }

} // namespace elemental
