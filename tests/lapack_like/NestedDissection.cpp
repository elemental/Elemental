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

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const Int n = Input("--n","size of n x n x n grid",30);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print graph?",false);
        const bool display = Input("--display","display graph?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        const Int numVertices = n*n*n;
        DistGraph graph( numVertices, comm );

        const Int firstLocalSource = graph.FirstLocalSource();
        const Int numLocalSources = graph.NumLocalSources();

        // Fill our portion of the graph of a 3D n x n x n 7-point stencil
        // in natural ordering: (x,y,z) at x + y*n + z*n*n
        if( commRank == 0 )
            Output("Filling local portion of graph");
        graph.Reserve( 7*numLocalSources );
        for( int iLocal=0; iLocal<numLocalSources; ++iLocal )
        {
            const int i = firstLocalSource + iLocal;
            const int x = i % n;
            const int y = (i/n) % n;
            const int z = i/(n*n);

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

        if( commRank == 0 )
            Output("Running nested dissection");
        ldl::DistNodeInfo info;
        ldl::DistSeparator sep;
        DistMap map;
        ldl::NestedDissection( graph, map, sep, info, ctrl );

        const int rootSepSize = info.size;
        // TODO: Print more than just the root separator size
        if( commRank == 0 )
            Output(rootSepSize," vertices in root separator");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
