/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalEnvironment.h"
#include "ElementalLAPACKInternal.h"
using namespace std;
using namespace Elemental;

static bool elementalInitializedMPI;
static stack<int> blocksizeStack;

void
Elemental::Init
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
    LAPACK::Internal::CreatePivotOp<float>();
    LAPACK::Internal::CreatePivotOp<double>();
#ifndef WITHOUT_COMPLEX
    LAPACK::Internal::CreatePivotOp<scomplex>();
    LAPACK::Internal::CreatePivotOp<dcomplex>();
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
Elemental::Finalize()
{
#ifndef RELEASE
    PushCallStack("Finalize");
#endif
    int finalized;
    MPI_Finalized( &finalized );
    if( finalized != 0 )
        cerr << "Warning: MPI was finalized before Elemental." << endl;

    LAPACK::Internal::DestroyPivotOp<float>();
    LAPACK::Internal::DestroyPivotOp<double>();
#ifndef WITHOUT_COMPLEX
    LAPACK::Internal::DestroyPivotOp<scomplex>();
    LAPACK::Internal::DestroyPivotOp<dcomplex>();
#endif
    if( ::elementalInitializedMPI && finalized == 0 )
        MPI_Finalize();
#ifndef RELEASE
    PopCallStack();
#endif
}

int 
Elemental::Blocksize()
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( 192 );
    return ::blocksizeStack.top(); 
}

void
Elemental::SetBlocksize( int blocksize )
{ 
    if( ::blocksizeStack.size() == 0 )
        ::blocksizeStack.push( blocksize );
    else
        ::blocksizeStack.top() = blocksize; 
}

void
Elemental::PushBlocksizeStack( int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void
Elemental::PopBlocksizeStack()
{ ::blocksizeStack.pop(); }

// If we are not in RELEASE mode, then implement wrappers for a CallStack
#ifndef RELEASE
static stack<string> callStack;

void
Elemental::PushCallStack( string name )
{ ::callStack.push( name ); }

void
Elemental::PopCallStack()
{ ::callStack.pop(); }

void
Elemental::DumpCallStack()
{
    while( ! ::callStack.empty() )
    {
        cerr << "[" << ::callStack.size() << "]: " << ::callStack.top() << endl;
        ::callStack.pop();
    }
}
#endif

