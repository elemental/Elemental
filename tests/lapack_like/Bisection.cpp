/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );

    try
    {
        const Int n = Input("--n","size of n x n x n grid",30);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const Int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const Int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const bool print = Input("--print","print graph?",false);
        const bool display = Input("--display","display graph?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;

        const Int numVertices = n*n*n;
        DistGraph graph( numVertices, comm );

        // Fill our portion of the graph of a 3D n x n x n 7-point stencil
        // in natural ordering: (x,y,z) at x + y*n + z*n*n
        const Int firstLocalSource = graph.FirstLocalSource();
        const Int numLocalSources = graph.NumLocalSources();
        graph.Reserve( 7*numLocalSources );
        for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
        {
            const Int i = firstLocalSource + iLocal;
            const Int x = i % n;
            const Int y = (i/n) % n;
            const Int z = i/(n*n);

            graph.QueueLocalConnection( iLocal, i );
            if( x != 0 )
                graph.QueueLocalConnection( iLocal, i-1 );
            if( x != n-1 )
                graph.QueueLocalConnection( iLocal, i+1 );
            if( y != 0 )
                graph.QueueLocalConnection( iLocal, i-n );
            if( y != n-1 )
                graph.QueueLocalConnection( iLocal, i+n );
            if( z != 0 )
                graph.QueueLocalConnection( iLocal, i-n*n );
            if( z != n-1 )
                graph.QueueLocalConnection( iLocal, i+n*n );
        }
        graph.ProcessQueues();
        if( display )
            Display( graph );
        if( print )
            Print( graph );

        if( commSize > 1 )
        {
            DistGraph child;
            DistMap map;
            bool haveLeftChild;
            const Int sepSize = 
                Bisect( graph, child, map, haveLeftChild, ctrl );

            int leftChildSize, rightChildSize;
            if( haveLeftChild )
            {
                leftChildSize = child.NumSources();
                rightChildSize = numVertices - leftChildSize - sepSize;
            }
            else
            {
                rightChildSize = child.NumSources();
                leftChildSize = numVertices - rightChildSize - sepSize;
            }
            if( commRank == 0 )
            {
                if( haveLeftChild )
                    Output
                    ("Root is on left with sizes: ",leftChildSize,",",
                     rightChildSize,",",sepSize);
                else
                    Output
                    ("Root is on right with sizes: ",leftChildSize,",",
                     rightChildSize,",",sepSize);
            }
        }
        else
        {
            // Turn the single-process DistGraph into a Graph
            Graph seqGraph( graph );

            Graph leftChild, rightChild;
            std::vector<Int> map;
            const Int sepSize = 
                Bisect( seqGraph, leftChild, rightChild, map, ctrl );

            const Int leftChildSize = leftChild.NumSources();
            const Int rightChildSize = rightChild.NumSources();
            Output
            ("Partition sizes were: ",leftChildSize,",",rightChildSize,",",
             sepSize);
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
