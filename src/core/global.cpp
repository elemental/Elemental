/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
#ifdef HAVE_OPENMP
        const int provided = mpi::QueryThread();
        if( provided != mpi::THREAD_MULTIPLE )
        {
            throw std::runtime_error
            ("MPI initialized with inadequate thread support for Elemental");
        }
#endif
        ::elemInitializedMpi = false;
    }

    // Queue a default algorithmic blocksize
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::blocksizeStack.push( 128 );

    // Build the default grid
    defaultGrid = new Grid( mpi::COMM_WORLD );

    // Build the pivot operations needed by the distributed LU
    CreatePivotOp<float>();
    CreatePivotOp<double>();
    CreatePivotOp<Complex<float> >();
    CreatePivotOp<Complex<double> >();

    // Seed the random number generators using Katzgrabber's approach
    // from "Random Numbers in Scientific Computing: An Introduction"
    // NOTE: srand no longer needed after C++11
    const unsigned rank = mpi::CommRank( mpi::COMM_WORLD );
    const long secs = time(NULL);
    const long seed = abs(((secs*181)*((rank-83)*359))%104729);
#ifdef WIN32
    srand( seed );
#else
    srand48( seed );
#endif
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
            DestroyPivotOp<float>();
            DestroyPivotOp<double>();
            DestroyPivotOp<Complex<float> >();
            DestroyPivotOp<Complex<double> >();

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
{ 
#ifdef HAVE_OPENMP
    if( omp_get_thread_num() != 0 )
        return;
#endif // HAVE_OPENMP
    ::callStack.push(s); 
}

void PopCallStack()
{ 
#ifdef HAVE_OPENMP
    if( omp_get_thread_num() != 0 )
        return;
#endif // HAVE_OPENMP
    ::callStack.pop(); 
}

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
