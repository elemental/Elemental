/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/environment.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

static bool elementalInitializedMPI;
static stack<int> blocksizeStack;

void
elemental::Init
( int* argc, char** argv[] )
{
#ifndef RELEASE
    PushCallStack("init");
#endif
    int initialized;
    MPI_Initialized( &initialized );
    if( initialized == 0 )
    {
        int finalized; 
        MPI_Finalized( &finalized );
        if( finalized != 0 )
        {
            cerr << "Cannot initialize Elemental after MPI_Finalize." << endl;
#ifndef RELEASE
            DumpCallStack();
#endif
            throw exception();
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
    srand( rank );
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
static stack<string> callStack;

void
elemental::PushCallStack( string s )
{ ::callStack.push(s); }

void
elemental::PopCallStack()
{ ::callStack.pop(); }

void
elemental::DumpCallStack()
{
    while( ! ::callStack.empty() )
    {
        cerr << "[" << ::callStack.size() << "]: " << ::callStack.top() << endl;
        ::callStack.pop();
    }
}
#endif

