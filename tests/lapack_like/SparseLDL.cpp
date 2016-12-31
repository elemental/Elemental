/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.

   Copyright (c) 2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

// TODO(poulson): Modernize this test driver

template<typename Field>
void TestSparseDirect
( Int n1,
  Int n2,
  Int n3,
  Int numRHS,
  bool solve2d,
  bool selInv,
  bool intraPiv,
  Int nbFact,
  Int nbSolve,
  bool natural,
  bool unpack,
  bool print,
  bool display,
  const BisectCtrl& ctrl,
  const El::Grid& grid )
{
    typedef Base<Field> Real;
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());

    const int N = n1*n2*n3;
    DistSparseMatrix<Field> A(grid);
    Laplacian( A, n1, n2, n3 );
    A *= -Field(1);
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

    OutputFromRoot(grid.Comm(),"Generating random X and forming Y := A X...");
    timer.Start();
    DistMultiVec<Field> X( N, numRHS, grid ), Y( N, numRHS, grid );
    MakeUniform( X );
    Zero( Y );
    Multiply( NORMAL, Field(1), A, X, Field(0), Y );
    Matrix<Real> YOrigNorms;
    ColumnTwoNorms( Y, YOrigNorms );
    mpi::Barrier( grid.Comm() );
    timer.Stop();
    OutputFromRoot(grid.Comm(),timer.Partial()," seconds");
    if( display )
    {
        Display( X, "X" );
        Display( Y, "Y" );
    }
    if( print )
    {
        Print( X, "X" );
        Print( Y, "Y" );
    }

    OutputFromRoot(grid.Comm(),"Running analysis...");
    timer.Start();
    const bool hermitian = true;
    DistSparseLDLFactorization<Field> sparseLDLFact;
    if( natural )
        sparseLDLFact.Initialize3DGridGraph
        ( n1, n2, n3, A, hermitian, ctrl );
    else
        sparseLDLFact.Initialize( A, hermitian, ctrl );
    mpi::Barrier( grid.Comm() );
    timer.Stop();
    OutputFromRoot(grid.Comm(),timer.Partial()," seconds");

    const Int rootSepSize = sparseLDLFact.NodeInfo().size;
    OutputFromRoot(grid.Comm(),rootSepSize," vertices in root separator\n");
    // TODO(poulson): Update the following for the new data structures
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

    // Unpack the ldl::DistFront into a sparse matrix
    if( unpack )
    {
        DistSparseMatrix<Field> APerm;
        sparseLDLFact.Front().Unpack
        ( APerm, sparseLDLFact.Separator(), sparseLDLFact.NodeInfo() );
        MakeSymmetric( LOWER, APerm );
        if( print )
            Print( APerm, "APerm" );
        if( display )
            Display( APerm, "APerm" );
    }

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

    OutputFromRoot(grid.Comm(),"Solving against Y...");
    SetBlocksize( nbSolve );
    mpi::Barrier( grid.Comm() );
    timer.Start();
    sparseLDLFact.Solve( Y );
    mpi::Barrier( grid.Comm() );
    const double solveTime = timer.Stop();
    const double localSolveGFlops = sparseLDLFact.LocalSolveGFlops( numRHS );
    const double solveGFlops = mpi::AllReduce( localSolveGFlops, grid.Comm() );
    const double solveSpeed = solveGFlops / factTime;
    OutputFromRoot(grid.Comm(),solveTime," seconds (",solveSpeed," GFlop/s)");

    OutputFromRoot(grid.Comm(),"Checking error in computed solution...");
    Matrix<Real> XNorms, YNorms;
    ColumnTwoNorms( X, XNorms );
    ColumnTwoNorms( Y, YNorms );
    Y -= X;
    Matrix<Real> errorNorms;
    ColumnTwoNorms( Y, errorNorms );
    for( int j=0; j<numRHS; ++j )
        OutputFromRoot
        (grid.Comm(),
         "Right-hand side ",j,":\n",Indent(),
         "|| x     ||_2 = ",XNorms.Get(j,0),"\n",Indent(),
         "|| error ||_2 = ",errorNorms.Get(j,0),"\n",Indent(),
         "|| A x   ||_2 = ",YOrigNorms.Get(j,0),"\n");
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
        const Int numRHS = Input("--numRHS","number of right-hand sides",5);
        const bool solve2d = Input("--solve2d","use 2d solve?",false);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool intraPiv = Input("--intraPiv","pivot within fronts?",false);
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
        const Int nbSolve = Input("--nbSolve","solve blocksize",96);
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool unpack = Input("--unpack","unpack frontal matrix?",true);
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
        const El::Grid grid(comm);

        // TODO(poulson): Call complex variants as well

        TestSparseDirect<float>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
        TestSparseDirect<double>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
#ifdef EL_HAVE_QD
        TestSparseDirect<DoubleDouble>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
        TestSparseDirect<QuadDouble>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
#endif
#ifdef EL_HAVE_QUAD
        TestSparseDirect<Quad>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSparseDirect<BigFloat>
        ( n1, n2, n3, numRHS, solve2d, selInv, intraPiv, nbFact, nbSolve,
          natural, unpack, print, display, ctrl, grid );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
