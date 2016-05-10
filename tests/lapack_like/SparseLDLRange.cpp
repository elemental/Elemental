/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",30);
        const Int n3 = Input("--n3","third grid dimension",30);
        const Int numRHSBeg = Input("--numRHSBeg","min number of rhs's",100);
        const Int numRHSInc = Input("--numRHSInc","stepsize for rhs's",100);
        const Int numRHSEnd = Input("--numRHSEnd","max number of rhs's",1000);
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
        const Int nbFact = Input("--nbFact","factorization blocksize",96);
        const Int nbSolveBeg = Input("--nbSolveBeg","min solve blocksize",96);
        const Int nbSolveInc = Input("--nbSolveInc","stepsize for bsize",16);
        const Int nbSolveEnd = Input("--nbSolveEnd","max solve blocksize",256);
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        const Int N = n1*n2*n3;
        DistSparseMatrix<C> A(comm);
        Laplacian( A, n1, n2, n3 );
        A *= -1;
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

        Timer timer;

        OutputFromRoot(comm,"Running nested dissection...");
        timer.Start();
        const auto& graph = A.DistGraph();
        ldl::DistNodeInfo info;
        ldl::DistSeparator sep;
        DistMap map, invMap;
        if( natural )
            ldl::NaturalNestedDissection
            ( n1, n2, n3, graph, map, sep, info, cutoff );
        else
            ldl::NestedDissection( graph, map, sep, info, ctrl );
        InvertMap( map, invMap );
        mpi::Barrier( comm );
        timer.Stop();
        OutputFromRoot(comm,timer.Partial()," seconds");

        const Int rootSepSize = info.size;
        OutputFromRoot(comm,rootSepSize," vertices in root separator\n");
        /*
        if( display )
        {
            ostringstream osBefore, osAfter;
            osBefore << "Structure before fact. on process " << commRank;
            DisplayLocal( info, false, osBefore.str() );
            osAfter << "Structure after fact. on process " << commRank;
            DisplayLocal( info, true, osAfter.str() );
        }
        */

        OutputFromRoot(comm,"Building DistSymmFront tree...");
        mpi::Barrier( comm );
        timer.Start();
        ldl::DistFront<C> front( A, map, sep, info, false );
        mpi::Barrier( comm );
        timer.Stop();
        OutputFromRoot(comm,timer.Partial()," seconds");

        // Memory usage before factorization
        const Int localEntriesBefore = front.NumLocalEntries();
        const Int minLocalEntriesBefore =
          mpi::AllReduce( localEntriesBefore, mpi::MIN, comm );
        const Int maxLocalEntriesBefore =
          mpi::AllReduce( localEntriesBefore, mpi::MAX, comm );
        const Int entriesBefore =
          mpi::AllReduce( localEntriesBefore, mpi::SUM, comm );
        OutputFromRoot
        (comm,
         "Memory usage before factorization: \n",Indent(),
         "  min entries:   ",minLocalEntriesBefore,"\n",Indent(),
         "  max entries:   ",maxLocalEntriesBefore,"\n",Indent(),
         "  total entries: ",entriesBefore,"\n");

        OutputFromRoot(comm,"Running LDL^T and redistribution...");
        SetBlocksize( nbFact );
        mpi::Barrier( comm );
        timer.Start();
        LDLFrontType type;
        if( solve2d )
        {
            if( intraPiv )
                type = ( selInv ? LDL_INTRAPIV_SELINV_2D : LDL_INTRAPIV_2D );
            else
                type = ( selInv ? LDL_SELINV_2D : LDL_2D );
        }
        else
        {
            if( intraPiv )
                type = ( selInv ? LDL_INTRAPIV_SELINV_1D : LDL_INTRAPIV_1D );
            else
                type = ( selInv ? LDL_SELINV_1D : LDL_1D );
        }
        LDL( info, front, type );
        mpi::Barrier( comm );
        const double factTime = timer.Stop();
        const double localFactGFlops = front.LocalFactorGFlops( selInv );
        const double factGFlops = mpi::AllReduce( localFactGFlops, comm );
        const double factSpeed = factGFlops / factTime;
        OutputFromRoot(comm,factTime," seconds, ",factSpeed," GFlop/s");

        // Memory usage after factorization
        const Int localEntriesAfter = front.NumLocalEntries();
        const Int minLocalEntriesAfter =
          mpi::AllReduce( localEntriesAfter, mpi::MIN, comm );
        const Int maxLocalEntriesAfter =
          mpi::AllReduce( localEntriesAfter, mpi::MAX, comm );
        const Int entriesAfter =
          mpi::AllReduce( localEntriesAfter, mpi::SUM, comm );
        OutputFromRoot
        (comm,
         "Memory usage after factorization: \n",Indent(),
         "  min entries:   ",minLocalEntriesAfter,"\n",Indent(),
         "  max entries:   ",maxLocalEntriesAfter,"\n",Indent(),
         "  total entries: ",entriesAfter,"\n");

        for( Int numRHS=numRHSBeg; numRHS<=numRHSEnd; numRHS+=numRHSInc )
        {
            const double localSolveGFlops = front.LocalSolveGFlops(numRHS);
            const double solveGFlops = mpi::AllReduce( localSolveGFlops, comm );

            DistMultiVec<C> Y( N, numRHS, comm );
            for( Int nbSolve=nbSolveBeg; nbSolve<=nbSolveEnd; 
                 nbSolve+=nbSolveInc )
            {
                MakeUniform( Y );
                SetBlocksize( nbSolve );
                OutputFromRoot(comm,"nbSolve=",nbSolve,"...");
                mpi::Barrier( comm );
                timer.Start();
                ldl::SolveAfter( invMap, info, front, Y );
                mpi::Barrier( comm );
                const double solveTime = timer.Stop();
                const double solveSpeed = solveGFlops / solveTime;
                OutputFromRoot
                (comm,solveTime," seconds (",solveSpeed," GFlop/s)");
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
