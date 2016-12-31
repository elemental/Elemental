/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.

   Copyright (c) 2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Field>
void TestSparseDirect
( Int n1,
  Int n2,
  Int n3,
  Int numRHSBeg,
  Int numRHSInc,
  Int numRHSEnd,
  bool intraPiv,
  bool solve2d,
  bool selInv,
  Int nbFact,
  Int nbSolveBeg,
  Int nbSolveInc,
  Int nbSolveEnd,
  bool natural,
  bool print,
  bool display,
  const BisectCtrl& ctrl,
  const Grid& grid )
{
    const Int N = n1*n2*n3;
    DistSparseMatrix<Field> A(grid);
    Laplacian( A, n1, n2, n3 );
    A *= Field(-1);
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
    const bool hermitian = true;

    Timer timer;

    DistSparseLDLFactorization<Field> sparseLDLFact;
    OutputFromRoot(grid.Comm(),"Running nested dissection...");
    timer.Start();
    if( natural )
        sparseLDLFact.Initialize3DGridGraph( n1, n2, n3, A, hermitian, ctrl );
    else
        sparseLDLFact.Initialize( A, hermitian, ctrl );
    mpi::Barrier( grid.Comm() );
    timer.Stop();
    OutputFromRoot(grid.Comm(),timer.Partial()," seconds");

    const Int rootSepSize = sparseLDLFact.NodeInfo().size;
    OutputFromRoot(grid.Comm(),rootSepSize," vertices in root separator\n");
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

    // Memory usage before factorization
    const Int localEntriesBefore = sparseLDLFact.NumLocalEntries();
    const Int minLocalEntriesBefore =
      mpi::AllReduce( localEntriesBefore, mpi::MIN, grid.Comm() );
    const Int maxLocalEntriesBefore =
      mpi::AllReduce( localEntriesBefore, mpi::MAX, grid.Comm() );
    const Int entriesBefore =
      mpi::AllReduce( localEntriesBefore, mpi::SUM, grid.Comm() );
    OutputFromRoot
    (grid.Comm(),
     "Memory usage before factorization: \n",Indent(),
     "  min entries:   ",minLocalEntriesBefore,"\n",Indent(),
     "  max entries:   ",maxLocalEntriesBefore,"\n",Indent(),
     "  total entries: ",entriesBefore,"\n");

    OutputFromRoot(grid.Comm(),"Running LDL^T and redistribution...");
    SetBlocksize( nbFact );
    mpi::Barrier( grid.Comm() );
    timer.Start();
    LDLFrontType type;
    if( solve2d )
    {
        if( intraPiv )
            type = selInv ? LDL_INTRAPIV_SELINV_2D : LDL_INTRAPIV_2D;
        else
            type = selInv ? LDL_SELINV_2D : LDL_2D;
    }
    else
    {
        if( intraPiv )
            type = selInv ? LDL_INTRAPIV_SELINV_1D : LDL_INTRAPIV_1D;
        else
            type = selInv ? LDL_SELINV_1D : LDL_1D;
    }
    sparseLDLFact.Factor( type );
    mpi::Barrier( grid.Comm() );
    const double factTime = timer.Stop();
    const double localFactGFlops = sparseLDLFact.LocalFactorGFlops( selInv );
    const double factGFlops = mpi::AllReduce( localFactGFlops, grid.Comm() );
    const double factSpeed = factGFlops / factTime;
    OutputFromRoot(grid.Comm(),factTime," seconds, ",factSpeed," GFlop/s");

    // Memory usage after factorization
    const Int localEntriesAfter = sparseLDLFact.NumLocalEntries();
    const Int minLocalEntriesAfter =
      mpi::AllReduce( localEntriesAfter, mpi::MIN, grid.Comm() );
    const Int maxLocalEntriesAfter =
      mpi::AllReduce( localEntriesAfter, mpi::MAX, grid.Comm() );
    const Int entriesAfter =
      mpi::AllReduce( localEntriesAfter, mpi::SUM, grid.Comm() );
    OutputFromRoot
    (grid.Comm(),
     "Memory usage after factorization: \n",Indent(),
     "  min entries:   ",minLocalEntriesAfter,"\n",Indent(),
     "  max entries:   ",maxLocalEntriesAfter,"\n",Indent(),
     "  total entries: ",entriesAfter,"\n");

    for( Int numRHS=numRHSBeg; numRHS<=numRHSEnd; numRHS+=numRHSInc )
    {
        const double localSolveGFlops = sparseLDLFact.LocalSolveGFlops(numRHS);
        const double solveGFlops =
          mpi::AllReduce( localSolveGFlops, grid.Comm() );

        DistMultiVec<Field> Y( N, numRHS, grid );
        for( Int nbSolve=nbSolveBeg; nbSolve<=nbSolveEnd;
             nbSolve+=nbSolveInc )
        {
            MakeUniform( Y );
            SetBlocksize( nbSolve );
            OutputFromRoot(grid.Comm(),"nbSolve=",nbSolve,"...");
            mpi::Barrier( grid.Comm() );
            timer.Start();
            sparseLDLFact.Solve( Y );
            mpi::Barrier( grid.Comm() );
            const double solveTime = timer.Stop();
            const double solveSpeed = solveGFlops / solveTime;
            OutputFromRoot
            (grid.Comm(),solveTime," seconds (",solveSpeed," GFlop/s)");
        }
    }
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int n1 = Input("--n1","first grid dimension",20);
        const Int n2 = Input("--n2","second grid dimension",15);
        const Int n3 = Input("--n3","third grid dimension",10);
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
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        const Grid grid( comm );

        TestSparseDirect<float>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<float>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );

        TestSparseDirect<double>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<double>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );

#ifdef EL_HAVE_QD
        TestSparseDirect<DoubleDouble>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<DoubleDouble>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );

        TestSparseDirect<QuadDouble>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<QuadDouble>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
#endif

#ifdef EL_HAVE_QUAD
        TestSparseDirect<Quad>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<Quad>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
#endif

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSparseDirect<BigFloat>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
        TestSparseDirect<Complex<BigFloat>>
        ( n1, n2, n3, numRHSBeg, numRHSInc, numRHSEnd, intraPiv, solve2d,
          selInv, nbFact, nbSolveBeg, nbSolveInc, nbSolveEnd, natural,
          print, display, ctrl, grid );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
