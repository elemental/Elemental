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

template<typename F>
void TestSparseDirect
( Int n1,
  Int n2,
  Int n3,
  Int numRHS,
  bool tryLDL,
  bool print,  
  bool display,
  const BisectCtrl& ctrl,
  mpi::Comm& comm )
{
    typedef Base<F> Real;
    OutputFromRoot(comm,"Testing with ",TypeName<F>());

    const int N = n1*n2*n3;
    DistSparseMatrix<F> A(comm);
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

    OutputFromRoot(comm,"Generating random vector X and forming Y := A X");
    timer.Start();
    DistMultiVec<F> X( N, numRHS, comm ), Y( N, numRHS, comm );
    MakeUniform( X );
    Zero( Y );
    Multiply( NORMAL, F(1), A, X, F(0), Y );
    Matrix<Real> YOrigNorms;
    ColumnTwoNorms( Y, YOrigNorms );
    mpi::Barrier( comm );
    timer.Stop();
    OutputFromRoot(comm,timer.Partial()," seconds");

    OutputFromRoot(comm,"Solving...");
    timer.Start();
    SymmetricSolve( A, Y, conjugate, tryLDL, ctrl );
    timer.Stop();
    OutputFromRoot(comm,timer.Partial()," seconds");

    OutputFromRoot(comm,"Checking error in computed solution...");
    Matrix<Real> XNorms, YNorms;
    ColumnTwoNorms( X, XNorms );
    ColumnTwoNorms( Y, YNorms );
    Y -= X;
    Matrix<Real> errorNorms;
    ColumnTwoNorms( Y, errorNorms );
    for( int j=0; j<numRHS; ++j )
        OutputFromRoot
        (comm,
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

        TestSparseDirect<float>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );
        TestSparseDirect<double>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );

#ifdef EL_HAVE_QD
        TestSparseDirect<DoubleDouble>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );
        TestSparseDirect<QuadDouble>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );
#endif

#ifdef EL_HAVE_QUAD
        TestSparseDirect<Quad>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );
#endif

#ifdef EL_HAVE_MPC
        TestSparseDirect<BigFloat>
        ( n1, n2, n3, numRHS, tryLDL, print, display, ctrl, comm );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
