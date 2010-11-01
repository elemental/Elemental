/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

namespace {
bool elementalInitializedMPI;
stack<int> blocksizeStack;
}

void
elemental::Init
( int* argc, char** argv[] )
{
#ifndef RELEASE
    PushCallStack("Init");
#endif
    int initialized;
    MPI_Initialized( &initialized );
    if( initialized == 0 )
    {
        int finalized; 
        MPI_Finalized( &finalized );
        if( finalized != 0 )
        {
            throw logic_error
            ( "Cannot initialize elemental after MPI_Finalize." );
        }
        MPI_Init( argc, argv );
        ::elementalInitializedMPI = true;
    }
    else
    {
        ::elementalInitializedMPI = false;
    }
    lapack::internal::CreatePivotOp<float>();
    lapack::internal::CreatePivotOp<double>();
#ifndef WITHOUT_COMPLEX
    lapack::internal::CreatePivotOp<scomplex>();
    lapack::internal::CreatePivotOp<dcomplex>();
#endif

    // Seed the random number generator with out rank
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    srand( rank+time(0) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    int finalized;
    MPI_Finalized( &finalized );
    if( finalized != 0 )
        cerr << "Warning: MPI was finalized before Elemental." << endl;

    lapack::internal::DestroyPivotOp<float>();
    lapack::internal::DestroyPivotOp<double>();
#ifndef WITHOUT_COMPLEX
    lapack::internal::DestroyPivotOp<scomplex>();
    lapack::internal::DestroyPivotOp<dcomplex>();
#endif
    if( ::elementalInitializedMPI && finalized == 0 )
        MPI_Finalize();
#ifndef RELEASE
    PopCallStack();
#endif
}

int 
elemental::Blocksize()
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( 192 );
    return ::blocksizeStack.top(); 
}

void
elemental::SetBlocksize( int blocksize )
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( blocksize );
    else
        ::blocksizeStack.top() = blocksize; 
}

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

