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
#include "elemental/environment.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

namespace {
bool elementalInitializedMpi;
stack<int> blocksizeStack;
}

void
elemental::Init
( int& argc, char**& argv )
{
    if( !mpi::Initialized() )
    {
        if( mpi::Finalized() )
        {
            throw logic_error
            ("Cannot initialize elemental after finalizing MPI.");
        }
#ifdef _OPENMP
        const int provided = 
            mpi::InitThread( argc, argv, mpi::THREAD_MULTIPLE );
        if( provided != mpi::THREAD_MULTIPLE )
        {
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
        }
#else
        mpi::Init( argc, argv );
#endif
        ::elementalInitializedMpi = true;
    }
    else
    {
        ::elementalInitializedMpi = false;
    }

    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::blocksizeStack.push( 128 );

    advanced::internal::CreatePivotOp<float>();
    advanced::internal::CreatePivotOp<double>();
#ifndef WITHOUT_COMPLEX
    advanced::internal::CreatePivotOp<scomplex>();
    advanced::internal::CreatePivotOp<dcomplex>();
#endif

    // Seed the parallel LCG
    plcg::UInt64 seed;
    seed.d[0] = time(0);
    seed.d[1] = time(0);
    const unsigned rank = mpi::CommRank( mpi::COMM_WORLD );
    const unsigned size = mpi::CommSize( mpi::COMM_WORLD );
    mpi::Broadcast( (char*)seed.d, 2*sizeof(unsigned), 0, mpi::COMM_WORLD );
    plcg::SeedParallelLcg( rank, size, seed );
}

void
elemental::Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    if( mpi::Finalized() )
        cerr << "Warning: MPI was finalized before Elemental." << endl;
    else if( ::elementalInitializedMpi )
    {
        advanced::internal::DestroyPivotOp<float>();
        advanced::internal::DestroyPivotOp<double>();
#ifndef WITHOUT_COMPLEX
        advanced::internal::DestroyPivotOp<scomplex>();
        advanced::internal::DestroyPivotOp<dcomplex>();
#endif
        mpi::Finalize();
    }
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

// If we are not in RELEASE mode, then implement wrappers for a CallStack
#ifndef RELEASE
namespace { stack<string> callStack; }

void
elemental::PushCallStack( string s )
{ ::callStack.push(s); }

void
elemental::PopCallStack()
{ ::callStack.pop(); }

void
elemental::DumpCallStack()
{
    ostringstream msg;
    while( ! ::callStack.empty() )
    {
        msg << "[" << ::callStack.size() << "]: " << ::callStack.top() << "\n"; 
        ::callStack.pop();
    }
    cerr << msg.str() << endl;;
}
#endif

