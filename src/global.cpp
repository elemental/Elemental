/*
   Copyright (c) 2009-2011, Jack Poulson
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
using namespace elemental;

//----------------------------------------------------------------------------//
// Variables for core routines                                                //
//----------------------------------------------------------------------------//
namespace {
bool initializedElemental = false;
bool elementalInitializedMpi;
std::stack<int> blocksizeStack;
Grid* defaultGrid = 0;
}

bool
elemental::Initialized()
{ return ::initializedElemental; }

void
elemental::Initialize( int& argc, char**& argv )
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
    advanced::internal::CreatePivotOp<float>();
    advanced::internal::CreatePivotOp<double>();
#ifndef WITHOUT_COMPLEX
    advanced::internal::CreatePivotOp<scomplex>();
    advanced::internal::CreatePivotOp<dcomplex>();
#endif

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

void
elemental::Finalize()
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
        advanced::internal::DestroyPivotOp<float>();
        advanced::internal::DestroyPivotOp<double>();
#ifndef WITHOUT_COMPLEX
        advanced::internal::DestroyPivotOp<scomplex>();
        advanced::internal::DestroyPivotOp<dcomplex>();
#endif

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

int
elemental::Blocksize()
{ return ::blocksizeStack.top(); }

void
elemental::SetBlocksize( int blocksize )
{ ::blocksizeStack.top() = blocksize; }

void
elemental::PushBlocksizeStack( int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void
elemental::PopBlocksizeStack()
{ ::blocksizeStack.pop(); }

const Grid&
elemental::DefaultGrid()
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
namespace { std::stack<std::string> callStack; }

void
elemental::PushCallStack( std::string s )
{ ::callStack.push(s); }

void
elemental::PopCallStack()
{ ::callStack.pop(); }

void
elemental::DumpCallStack()
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

//----------------------------------------------------------------------------//
// Variables for basic routines                                               //
//----------------------------------------------------------------------------//
namespace {
int localHemvFloatBlocksize = 64;
int localHemvDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localHemvComplexFloatBlocksize = 64;
int localHemvComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX

int localSymvFloatBlocksize = 64;
int localSymvDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localSymvComplexFloatBlocksize = 64;
int localSymvComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX

int localTriangularRank2KFloatBlocksize = 64;
int localTriangularRank2KDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localTriangularRank2KComplexFloatBlocksize = 64;
int localTriangularRank2KComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX

int localTriangularRankKFloatBlocksize = 64;
int localTriangularRankKDoubleBlocksize = 64;
#ifndef WITHOUT_COMPLEX
int localTriangularRankKComplexFloatBlocksize = 64;
int localTriangularRankKComplexDoubleBlocksize = 64;
#endif // WITHOUT_COMPLEX
}

template<>
void
elemental::basic::SetLocalHemvBlocksize<float>
( int blocksize )
{ ::localHemvFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalHemvBlocksize<double>
( int blocksize )
{ ::localHemvDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::basic::SetLocalHemvBlocksize< std::complex<float> >
( int blocksize )
{ ::localHemvComplexFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalHemvBlocksize< std::complex<double> >
( int blocksize )
{ ::localHemvComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::basic::LocalHemvBlocksize<float>()
{ return ::localHemvFloatBlocksize; }

template<>
int
elemental::basic::LocalHemvBlocksize<double>()
{ return ::localHemvDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::basic::LocalHemvBlocksize<scomplex>()
{ return ::localHemvComplexFloatBlocksize; }

template<>
int
elemental::basic::LocalHemvBlocksize<dcomplex>()
{ return ::localHemvComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

template<>
void
elemental::basic::SetLocalSymvBlocksize<float>
( int blocksize )
{ ::localSymvFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalSymvBlocksize<double>
( int blocksize )
{ ::localSymvDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::basic::SetLocalSymvBlocksize< std::complex<float> >
( int blocksize )
{ ::localSymvComplexFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalSymvBlocksize< std::complex<double> >
( int blocksize )
{ ::localSymvComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::basic::LocalSymvBlocksize<float>()
{ return ::localSymvFloatBlocksize; }

template<>
int
elemental::basic::LocalSymvBlocksize<double>()
{ return ::localSymvDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::basic::LocalSymvBlocksize<scomplex>()
{ return ::localSymvComplexFloatBlocksize; }

template<>
int
elemental::basic::LocalSymvBlocksize<dcomplex>()
{ return ::localSymvComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

template<>
void
elemental::basic::SetLocalTriangularRank2KBlocksize<float>
( int blocksize )
{ ::localTriangularRank2KFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRank2KBlocksize<double>
( int blocksize )
{ ::localTriangularRank2KDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::basic::SetLocalTriangularRank2KBlocksize< std::complex<float> >
( int blocksize )
{ ::localTriangularRank2KComplexFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRank2KBlocksize< std::complex<double> >
( int blocksize )
{ ::localTriangularRank2KComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::basic::LocalTriangularRank2KBlocksize<float>()
{ return ::localTriangularRank2KFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRank2KBlocksize<double>()
{ return ::localTriangularRank2KDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::basic::LocalTriangularRank2KBlocksize<scomplex>()
{ return ::localTriangularRank2KComplexFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRank2KBlocksize<dcomplex>()
{ return ::localTriangularRank2KComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize<float>
( int blocksize )
{ ::localTriangularRankKFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize<double>
( int blocksize )
{ ::localTriangularRankKDoubleBlocksize = blocksize; }

#ifndef WITHOUT_COMPLEX
template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize< std::complex<float> >
( int blocksize )
{ ::localTriangularRankKComplexFloatBlocksize = blocksize; }

template<>
void
elemental::basic::SetLocalTriangularRankKBlocksize< std::complex<double> >
( int blocksize )
{ ::localTriangularRankKComplexDoubleBlocksize = blocksize; }
#endif // WITHOUT_COMPLEX

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<float>()
{ return ::localTriangularRankKFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<double>()
{ return ::localTriangularRankKDoubleBlocksize; }

#ifndef WITHOUT_COMPLEX
template<>
int
elemental::basic::LocalTriangularRankKBlocksize<scomplex>()
{ return ::localTriangularRankKComplexFloatBlocksize; }

template<>
int
elemental::basic::LocalTriangularRankKBlocksize<dcomplex>()
{ return ::localTriangularRankKComplexDoubleBlocksize; }
#endif // WITHOUT_COMPLEX

//----------------------------------------------------------------------------//
// Variables for advanced routines                                            //
//----------------------------------------------------------------------------//
namespace {
HermitianTridiagApproach tridiagApproach = HERMITIAN_TRIDIAG_DEFAULT;
GridOrder gridOrder = ROW_MAJOR;
}

void
elemental::advanced::SetHermitianTridiagApproach
( HermitianTridiagApproach approach )
{ ::tridiagApproach = approach; }

HermitianTridiagApproach
elemental::advanced::GetHermitianTridiagApproach()
{ return ::tridiagApproach; }

void
elemental::advanced::SetHermitianTridiagGridOrder( GridOrder order )
{ ::gridOrder = order; }

GridOrder
elemental::advanced::GetHermitianTridiagGridOrder()
{ return ::gridOrder; }

