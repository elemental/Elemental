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
        const int numRhsBeg = Input("--numRhsBeg","min number of rhs's",100);
        const int numRhsInc = Input("--numRhsInc","stepsize for rhs's",100);
        const int numRhsEnd = Input("--numRhsEnd","max number of rhs's",1000);
        const bool intraPiv = Input("--intraPiv","frontal pivoting?",false);
        const bool solve2d = Input("--solve2d","use 2d solve?",false);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool natural = Input("--natural","analytical nested-diss?",true);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const int nbFact = Input("--nbFact","factorization blocksize",96);
        const int nbSolveBeg = Input("--nbSolveBeg","min solve blocksize",96);
        const int nbSolveInc = Input("--nbSolveInc","stepsize for bsize",16);
        const int nbSolveEnd = Input("--nbSolveEnd","max solve blocksize",256);
        const int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        const int N = n1*n2*n3;
        DistSparseMatrix<Complex<double> > A( N, comm );

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

            A.QueueUpdate( i, i, 6. );
            if( x != 0 )
                A.QueueUpdate( i, i-1, -1. );
            if( x != n1-1 )
                A.QueueUpdate( i, i+1, -1. );
            if( y != 0 )
                A.QueueUpdate( i, i-n1, -1. );
            if( y != n2-1 )
                A.QueueUpdate( i, i+n1, -1. );
            if( z != 0 )
                A.QueueUpdate( i, i-n1*n2, -1. );
            if( z != n3-1 )
                A.QueueUpdate( i, i+n1*n2, -1. );
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
        if( natural )
        {
            NaturalNestedDissection
            ( n1, n2, n3, graph, map, sepTree, info, cutoff );
        }
        else
        {
            NestedDissection
            ( graph, map, sepTree, info, 
              sequential, numDistSeps, numSeqSeps, cutoff );
        }
        map.FormInverse( inverseMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << nestedStop-nestedStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            const int distNodes = info.distNodes.size();
            const int localNodes = info.localNodes.size();
            const int rootSepSize = info.distNodes.back().size;
            std::cout << "\n"
                      << "On the root process:\n"
                      << "-----------------------------------------\n"
                      << localNodes << " local nodes\n"
                      << distNodes  << " distributed nodes\n"
                      << rootSepSize << " vertices in root separator\n"
                      << std::endl;
        }
        if( display )
        {
            std::ostringstream osBefore, osAfter;
            osBefore << "Structure before fact. on process " << commRank;
            osAfter << "Structure after fact. on process " << commRank;
            DisplayLocal( info, false, osBefore.str() );
            DisplayLocal( info, true, osAfter.str() );
        }

        if( commRank == 0 )
        {
            std::cout << "Building DistSymmFrontTree...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFrontTree<Complex<double>> 
            frontTree( A, map, sepTree, info, false );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << buildStop-buildStart << " seconds"
                      << std::endl;

        double localEntries, minLocalEntries, maxLocalEntries, globalEntries;
        frontTree.MemoryInfo
        ( localEntries, minLocalEntries, maxLocalEntries, globalEntries );
        double localFactFlops, minLocalFactFlops, maxLocalFactFlops, 
               globalFactFlops;
        frontTree.FactorizationWork
        ( localFactFlops, minLocalFactFlops, maxLocalFactFlops, 
          globalFactFlops, selInv );
        if( commRank == 0 )
        {
            std::cout 
              << "Original memory usage for fronts...\n"
              << "  min local: " << minLocalEntries*2*sizeof(double)/1e6 
              << " MB\n"
              << "  max local: " << maxLocalEntries*2*sizeof(double)/1e6 
              << " MB\n"
              << "  global:    " << globalEntries*2*sizeof(double)/1e6
              << " MB\n"
              << "\n"
              << "Factorization (and possibly sel-inv) work...\n"
              << "  min local: " << minLocalFactFlops/1.e9 << " GFlops\n"
              << "  max local: " << maxLocalFactFlops/1.e9 << " GFlops\n"
              << "  global:    " << globalFactFlops/1.e9 << " GFlops\n"
              << std::endl;
        }

        if( commRank == 0 )
        {
            std::cout << "Running LDL^T and redistribution...";
            std::cout.flush();
        }
        SetBlocksize( nbFact );
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        if( solve2d )
        {
            if( intraPiv )
            {
                if( selInv )    
                    LDL( info, frontTree, LDL_INTRAPIV_SELINV_2D );
                else
                    LDL( info, frontTree, LDL_INTRAPIV_2D );
            }
            else
            {
                if( selInv )
                    LDL( info, frontTree, LDL_SELINV_2D );
                else
                    LDL( info, frontTree, LDL_2D );
            }
        }
        else
        {
            if( intraPiv )
            {
                if( selInv )
                    LDL( info, frontTree, LDL_INTRAPIV_SELINV_2D );
                else
                    LDL( info, frontTree, LDL_INTRAPIV_1D );
            }
            else
            {
                if( selInv )
                    LDL( info, frontTree, LDL_SELINV_2D );
                else
                    LDL( info, frontTree, LDL_1D );
            }
        }
        mpi::Barrier( comm );
        const double ldlStop = mpi::Time();
        const double factTime = ldlStop - ldlStart;
        const double factGFlops = globalFactFlops/(1.e9*factTime);
        if( commRank == 0 )
            std::cout << "done, " << factTime << " seconds, " 
                      << factGFlops << " GFlop/s" << std::endl;

        if( commRank == 0 )
            std::cout << "Memory usage for fronts after factorization..."
                      << std::endl;
        frontTree.MemoryInfo
        ( localEntries, minLocalEntries, maxLocalEntries, globalEntries );
        if( commRank == 0 )
        {
            std::cout << "  min local: " << minLocalEntries*2*sizeof(double)/1e6
                      << " MB\n"
                      << "  max local: " << maxLocalEntries*2*sizeof(double)/1e6
                      << " MB\n"
                      << "  global:    " << globalEntries*2*sizeof(double)/1e6
                      << " MB\n"
                      << std::endl;
        }

        for( int numRhs=numRhsBeg; numRhs<=numRhsEnd; numRhs+=numRhsInc )
        {
            double localSolveFlops, minLocalSolveFlops, maxLocalSolveFlops,
                   globalSolveFlops;
            frontTree.SolveWork
            ( localSolveFlops, minLocalSolveFlops, maxLocalSolveFlops,
              globalSolveFlops, numRhs );
            if( commRank == 0 )
            {
                std::cout
                  << "Solve with " << numRhs << " right-hand sides...\n"
                  << "  min local: " << minLocalSolveFlops/1.e9 << " GFlops\n"
                  << "  max local: " << maxLocalSolveFlops/1.e9 << " GFlops\n"
                  << "  global:    " << globalSolveFlops/1.e9 << " GFlops\n"
                  << std::endl;
            }

            DistMultiVec<Complex<double> > Y( N, numRhs, comm );
            for( int nbSolve=nbSolveBeg; nbSolve<=nbSolveEnd; 
                 nbSolve+=nbSolveInc )
            {
                MakeUniform( Y );
                SetBlocksize( nbSolve );
                if( commRank == 0 )
                {
                    std::cout << "  nbSolve=" << nbSolve << "...";
                    std::cout.flush();
                }
                double solveStart, solveStop;
                if( solve2d )
                {
                    DistNodalMatrix<Complex<double> > YNodal;
                    YNodal.Pull( inverseMap, info, Y );
                    mpi::Barrier( comm );
                    solveStart = mpi::Time();
                    Solve( info, frontTree, YNodal );
                    mpi::Barrier( comm );
                    solveStop = mpi::Time();
                    YNodal.Push( inverseMap, info, Y );
                }
                else
                {
                    DistNodalMultiVec<Complex<double> > YNodal;
                    YNodal.Pull( inverseMap, info, Y );
                    mpi::Barrier( comm );
                    solveStart = mpi::Time();
                    Solve( info, frontTree, YNodal );
                    mpi::Barrier( comm );
                    solveStop = mpi::Time();
                    YNodal.Push( inverseMap, info, Y );
                }
                const double solveTime = solveStop - solveStart;
                const double solveGFlops = globalSolveFlops/(1.e9*solveTime);
                if( commRank == 0 )
                    std::cout << "done, " << solveTime << " seconds, "
                              << solveGFlops << " GFlop/s" << std::endl;
            }
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
