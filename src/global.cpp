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
int numElemInits = 0;
bool elemInitializedMpi;
std::stack<int> blocksizeStack;
elem::Grid* defaultGrid = 0;
elem::MpiArgs* args = 0;

// Debugging
#ifndef RELEASE
std::stack<std::string> callStack;
#endif

// Tuning parameters for basic routines
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
using namespace elem;
HermitianTridiagApproach tridiagApproach = HERMITIAN_TRIDIAG_DEFAULT;
GridOrder gridOrder = ROW_MAJOR;
}

namespace elem {

bool Initialized()
{ return ::numElemInits > 0; }

void Initialize( int& argc, char**& argv )
{
    if( ::numElemInits > 0 )
    {
        ++::numElemInits;
        return;
    }

    ::args = new MpiArgs( argc, argv );

    ::numElemInits = 1;
    if( !mpi::Initialized() )
    {
        if( mpi::Finalized() )
        {
            throw std::logic_error
            ("Cannot initialize elemental after finalizing MPI");
        }
#ifdef HAVE_OPENMP
        const int provided = 
            mpi::InitializeThread
            ( argc, argv, mpi::THREAD_MULTIPLE );
        const int commRank = mpi::CommRank( mpi::COMM_WORLD );
        if( provided != mpi::THREAD_MULTIPLE && commRank == 0 )
        {
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
        }
#else
        mpi::Initialize( argc, argv );
#endif
        ::elemInitializedMpi = true;
    }
    else
    {
        ::elemInitializedMpi = false;
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
    internal::CreatePivotOp<Complex<float> >();
    internal::CreatePivotOp<Complex<double> >();

    // Seed the parallel random number generator, PLCG
    plcg::UInt64 seed;
    seed.d[0] = time(0);
    seed.d[1] = time(0);
    const unsigned rank = mpi::CommRank( mpi::COMM_WORLD );
    const unsigned size = mpi::CommSize( mpi::COMM_WORLD );
    mpi::Broadcast
    ( (byte*)seed.d, 2*sizeof(unsigned), 0, mpi::COMM_WORLD );
    plcg::SeedParallelLcg( rank, size, seed );
}

void Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    if( ::numElemInits <= 0 )
        throw std::logic_error("Finalized Elemental more than initialized");
    --::numElemInits;

    if( mpi::Finalized() )
    {
        std::cerr << "Warning: MPI was finalized before Elemental." 
                  << std::endl;
    }
    if( ::numElemInits == 0 )
    {
        delete ::args;
        ::args = 0;

        if( ::elemInitializedMpi )
        {
            // Destroy the pivot ops needed by the distributed LU
            internal::DestroyPivotOp<float>();
            internal::DestroyPivotOp<double>();
            internal::DestroyPivotOp<Complex<float> >();
            internal::DestroyPivotOp<Complex<double> >();

            // Delete the default grid
            delete ::defaultGrid;
            ::defaultGrid = 0;

            mpi::Finalize();
        }

        delete ::defaultGrid;
        ::defaultGrid = 0;
        while( ! ::blocksizeStack.empty() )
            ::blocksizeStack.pop();
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

MpiArgs& GetArgs()
{ 
    if( args == 0 )
        throw std::runtime_error("No available instance of MpiArgs");
    return *::args; 
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
void SetLocalHemvBlocksize<Complex<float> >( int blocksize )
{ ::localHemvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalHemvBlocksize<Complex<double> >( int blocksize )
{ ::localHemvComplexDoubleBlocksize = blocksize; }

template<>
int LocalHemvBlocksize<float>()
{ return ::localHemvFloatBlocksize; }

template<>
int LocalHemvBlocksize<double>()
{ return ::localHemvDoubleBlocksize; }

template<>
int LocalHemvBlocksize<Complex<float> >()
{ return ::localHemvComplexFloatBlocksize; }

template<>
int LocalHemvBlocksize<Complex<double> >()
{ return ::localHemvComplexDoubleBlocksize; }

template<>
void SetLocalSymvBlocksize<float>( int blocksize )
{ ::localSymvFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<double>( int blocksize )
{ ::localSymvDoubleBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<float> >( int blocksize )
{ ::localSymvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<double> >( int blocksize )
{ ::localSymvComplexDoubleBlocksize = blocksize; }

template<>
int LocalSymvBlocksize<float>()
{ return ::localSymvFloatBlocksize; }

template<>
int LocalSymvBlocksize<double>()
{ return ::localSymvDoubleBlocksize; }

template<>
int LocalSymvBlocksize<Complex<float> >()
{ return ::localSymvComplexFloatBlocksize; }

template<>
int LocalSymvBlocksize<Complex<double> >()
{ return ::localSymvComplexDoubleBlocksize; }

template<>
void SetLocalTrr2kBlocksize<float>( int blocksize )
{ ::localTrr2kFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<double>( int blocksize )
{ ::localTrr2kDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<float> >( int blocksize )
{ ::localTrr2kComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<double> >( int blocksize )
{ ::localTrr2kComplexDoubleBlocksize = blocksize; }

template<>
int LocalTrr2kBlocksize<float>()
{ return ::localTrr2kFloatBlocksize; }

template<>
int LocalTrr2kBlocksize<double>()
{ return ::localTrr2kDoubleBlocksize; }

template<>
int LocalTrr2kBlocksize<Complex<float> >()
{ return ::localTrr2kComplexFloatBlocksize; }

template<>
int LocalTrr2kBlocksize<Complex<double> >()
{ return ::localTrr2kComplexDoubleBlocksize; }

template<>
void SetLocalTrrkBlocksize<float>( int blocksize )
{ ::localTrrkFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<double>( int blocksize )
{ ::localTrrkDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<float> >( int blocksize )
{ ::localTrrkComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<double> >( int blocksize )
{ ::localTrrkComplexDoubleBlocksize = blocksize; }

template<>
int LocalTrrkBlocksize<float>()
{ return ::localTrrkFloatBlocksize; }

template<>
int LocalTrrkBlocksize<double>()
{ return ::localTrrkDoubleBlocksize; }

template<>
int LocalTrrkBlocksize<Complex<float> >()
{ return ::localTrrkComplexFloatBlocksize; }

template<>
int LocalTrrkBlocksize<Complex<double> >()
{ return ::localTrrkComplexDoubleBlocksize; }

void SetHermitianTridiagApproach( HermitianTridiagApproach approach )
{ ::tridiagApproach = approach; }

HermitianTridiagApproach GetHermitianTridiagApproach()
{ return ::tridiagApproach; }

void SetHermitianTridiagGridOrder( GridOrder order )
{ ::gridOrder = order; }

GridOrder GetHermitianTridiagGridOrder()
{ return ::gridOrder; }

} // namespace elem
