/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
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
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const int n1 = Input("--n1","first grid dimension",30);
        const int n2 = Input("--n2","second grid dimension",30);
        const int n3 = Input("--n3","third grid dimension",30);
        const int numRepeats = Input
            ("--numRepeats","number of repeated factorizations",5);
        const bool intraPiv = Input("--intraPiv","frontal pivoting?",false);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of partitions to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of partitions to try per sequential partition",1);
        const int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        const int N = n1*n2*n3;
        DistSparseMatrix<double> A( N, comm );

        // Fill our portion of the 3D negative Laplacian using a n1 x n2 x n3
        // 7-point stencil in natural ordering: (x,y,z) at x + y*n1 + z*n1*n2
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        const double fillStart = mpi::Time();
        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();
        A.Reserve( 7*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n1;
            const int y = (i/n1) % n2;
            const int z = i/(n1*n2);

            A.QueueLocalUpdate( iLocal, i, 6. );
            if( x != 0 )
                A.QueueLocalUpdate( iLocal, i-1, -1. );
            if( x != n1-1 )
                A.QueueLocalUpdate( iLocal, i+1, -1. );
            if( y != 0 )
                A.QueueLocalUpdate( iLocal, i-n1, -1. );
            if( y != n2-1 )
                A.QueueLocalUpdate( iLocal, i+n1, -1. );
            if( z != 0 )
                A.QueueLocalUpdate( iLocal, i-n1*n2, -1. );
            if( z != n3-1 )
                A.QueueLocalUpdate( iLocal, i+n1*n2, -1. );
        } 
        A.MakeConsistent();
        mpi::Barrier( comm );
        const double fillStop =  mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << fillStop-fillStart << " seconds" 
                      << std::endl;
        if( display )
        {
            Display( A );
            Display( A.DistGraph() );
        }
        if( print )
        {
            Print( A );
            Print( A.DistGraph() );
        }

        if( commRank == 0 )
        {
            std::cout << "Running nested dissection...";
            std::cout.flush();
        }
        const double nestedStart = mpi::Time();
        const DistGraph& graph = A.DistGraph();
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        DistMap map, inverseMap;
        NestedDissection
        ( graph, map, sepTree, info, 
          sequential, numDistSeps, numSeqSeps, cutoff );
        map.FormInverse( inverseMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << nestedStop-nestedStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            const int numDistNodes = info.distNodes.size();
            const int numLocalNodes = info.localNodes.size();
            const int rootSepSize = info.distNodes.back().size;
            std::cout << "\n"
                      << "On the root process:\n"
                      << "-----------------------------------------\n"
                      << numLocalNodes << " local nodes\n"
                      << numDistNodes  << " distributed nodes\n"
                      << rootSepSize << " vertices in root separator\n"
                      << std::endl;
        }

        if( commRank == 0 )
        {
            std::cout << "Building DistSymmFrontTree...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFrontTree<double> frontTree( A, map, sepTree, info, false );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << buildStop-buildStart << " seconds"
                      << std::endl;

        for( int repeat=0; repeat<numRepeats; ++repeat )
        {
            if( repeat != 0 )
            {
                // Reset to an unfactored, implicitly symmetric frontal tree
                if( commRank == 0 )
                    std::cout << "Resetting frontal tree." << std::endl;
                ChangeFrontType( frontTree, SYMM_2D );

                // Randomize the fronts
                if( commRank == 0 )
                    std::cout << "Randomizing fronts." << std::endl;
                const int numDistFronts = frontTree.distFronts.size();
                const int numLocalFronts = frontTree.localFronts.size();
                for( int s=0; s<numLocalFronts; ++s )
                    MakeUniform( frontTree.localFronts[s].frontL );
                for( int s=1; s<numDistFronts; ++s )
                    MakeUniform( frontTree.distFronts[s].front2dL );
            }

            if( commRank == 0 )
            {
                std::cout << "Running LDL^T and redistribution...";
                std::cout.flush();
            }
            mpi::Barrier( comm );
            const double ldlStart = mpi::Time();
            if( intraPiv )
                LDL( info, frontTree, LDL_INTRAPIV_1D );
            else
                LDL( info, frontTree, LDL_1D );
            mpi::Barrier( comm );
            const double ldlStop = mpi::Time();
            if( commRank == 0 )
                std::cout << "done, " << ldlStop-ldlStart << " seconds" 
                          << std::endl;

            if( commRank == 0 )
            {
                std::cout << "Solving against random right-hand side...";
                std::cout.flush();
            }
            const double solveStart = mpi::Time();
            DistMultiVec<double> y( N, 1, comm );
            MakeUniform( y );
            DistNodalMultiVec<double> yNodal;
            yNodal.Pull( inverseMap, info, y );
            Solve( info, frontTree, yNodal );
            yNodal.Push( inverseMap, info, y );
            mpi::Barrier( comm );
            const double solveStop = mpi::Time();
            if( commRank == 0 )
                std::cout << "done, " << solveStop-solveStart << " seconds"
                          << std::endl;

            // TODO: Check residual error
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
