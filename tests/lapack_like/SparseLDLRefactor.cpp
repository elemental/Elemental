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
void MakeFrontsUniform( ldl::Front<Field>& front )
{
    MakeUniform( front.LDense );
    for( const auto& child : front.children )
        MakeFrontsUniform( *child );
}

template<typename Field>
void MakeFrontsUniform( ldl::DistFront<Field>& front )
{
    MakeUniform( front.L2D );
    if( front.child.get() != nullptr )
        MakeFrontsUniform( *front.child );
    else
        MakeFrontsUniform( *front.duplicate );
}

template<typename Field>
void TestSparseDirect
( Int n1,
  Int n2,
  Int n3,
  Int numRepeats,
  bool intraPiv,
  bool print,
  bool display,
  const BisectCtrl& ctrl,
  const El::Grid& grid )
{
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

    const bool hermitian = false;
    DistSparseLDLFactorization<Field> sparseLDLFact;
    Timer timer;
    mpi::Barrier( grid.Comm() );
    OutputFromRoot(grid.Comm(),"Building ldl::DistFront tree...");
    timer.Start();
    sparseLDLFact.Initialize( A, hermitian, ctrl );
    mpi::Barrier( grid.Comm() );
    timer.Stop();
    OutputFromRoot(grid.Comm(),timer.Partial()," seconds");

    const Int rootSepSize = sparseLDLFact.NodeInfo().size;
    OutputFromRoot(grid.Comm(),rootSepSize," vertices in root separator\n");

    for( Int repeat=0; repeat<numRepeats; ++repeat )
    {
        if( repeat != 0 )
        {
            sparseLDLFact.ChangeFrontType( SYMM_2D );
            MakeFrontsUniform( sparseLDLFact.Front() );
        }

        OutputFromRoot(grid.Comm(),"Running LDL^T and redistribution...");
        mpi::Barrier( grid.Comm() );
        timer.Start();
        if( intraPiv )
            sparseLDLFact.Factor( LDL_INTRAPIV_1D );
        else
            sparseLDLFact.Factor( LDL_1D );
        mpi::Barrier( grid.Comm() );
        timer.Stop();
        OutputFromRoot(grid.Comm(),timer.Partial()," seconds");

        OutputFromRoot(grid.Comm(),"Solving against random right-hand side...");
        DistMultiVec<Field> y( N, 1, grid );
        MakeUniform( y );
        timer.Start();
        sparseLDLFact.Solve( y );
        mpi::Barrier( grid.Comm() );
        timer.Stop();
        OutputFromRoot(grid.Comm(),"Time = ",timer.Partial()," seconds");

        // TODO(poulson): Check residual error
    }
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",30);
        const Int n3 = Input("--n3","third grid dimension",30);
        const Int numRepeats = Input
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
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        const El::Grid grid( comm );

        TestSparseDirect<float>
        ( n1, n2, n3, numRepeats, intraPiv, print, display, ctrl, grid );
        TestSparseDirect<double>
        ( n1, n2, n3, numRepeats, intraPiv, print, display, ctrl, grid );

        // TODO(poulson): Test more datatypes? It makes the runtime much longer
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
