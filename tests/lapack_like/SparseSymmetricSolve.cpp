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
  Int numRHS,
  bool tryLDL,
  bool print,
  bool display,
  const BisectCtrl& ctrl,
  const Grid& grid )
{
    typedef Base<Field> Real;
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());

    const int N = n1*n2*n3;
    DistSparseMatrix<Field> A(grid);
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
    const bool conjugate = false;

    Timer timer;

    OutputFromRoot
    (grid.Comm(),"Generating random vector X and forming Y := A X");
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

    OutputFromRoot(grid.Comm(),"Solving...");
    timer.Start();
    SymmetricSolve( A, Y, conjugate, tryLDL, ctrl );
    timer.Stop();
    OutputFromRoot(grid.Comm(),timer.Partial()," seconds");

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
         "Right-hand side ",j,"\n",Indent(),
         "------------------------------------------\n",Indent(),
         "|| x     ||_2 = ",XNorms.Get(j,0),"\n",Indent(),
         "|| xComp ||_2 = ",YNorms.Get(j,0),"\n",Indent(),
         "|| A x   ||_2 = ",YOrigNorms.Get(j,0),"\n",Indent(),
         "|| error ||_2 = ",errorNorms.Get(j,0),"\n");
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",20);
        const Int n3 = Input("--n3","third grid dimension",10);
        const Int numRHS = Input("--numRHS","number of right-hand sides",5);
        const bool tryLDL = Input("--tryLDL","try LDL?",true);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        // TODO(poulson): Test complex variants?

        const Grid grid( comm );

        TestSparseDirect<float>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );
        TestSparseDirect<double>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );

#ifdef EL_HAVE_QD
        TestSparseDirect<DoubleDouble>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );
        TestSparseDirect<QuadDouble>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );
#endif

#ifdef EL_HAVE_QUAD
        TestSparseDirect<Quad>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );
#endif

#ifdef EL_HAVE_MPC
        TestSparseDirect<BigFloat>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, grid );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
