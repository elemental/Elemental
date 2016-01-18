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

// TODO: Extend this driver to test with several different datatypes

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",30);
        const Int n3 = Input("--n3","third grid dimension",30);
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

        const int N = n1*n2*n3;
        DistSparseMatrix<double> A(comm);
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

        if( commRank == 0 )
            Output("Generating random vector X and forming Y := A X");
        const double multiplyStart = mpi::Time();
        DistMultiVec<double> X( N, numRHS, comm ), Y( N, numRHS, comm );
        MakeUniform( X );
        Zero( Y );
        Multiply( NORMAL, 1., A, X, 0., Y );
        Matrix<double> YOrigNorms;
        ColumnTwoNorms( Y, YOrigNorms );
        mpi::Barrier( comm );
        const double multiplyStop = mpi::Time();
        if( commRank == 0 )
            Output(multiplyStop-multiplyStart," seconds");

        if( commRank == 0 )
            Output("Solving...");
        const double solveStart = mpi::Time();
        SymmetricSolve( A, Y, conjugate, tryLDL, ctrl );
        const double solveStop = mpi::Time();
        if( commRank == 0 )
            Output(solveStop-solveStart," seconds");

        if( commRank == 0 )
            Output("Checking error in computed solution...");
        Matrix<double> XNorms, YNorms;
        ColumnTwoNorms( X, XNorms );
        ColumnTwoNorms( Y, YNorms );
        Y -= X;
        Matrix<double> errorNorms;
        ColumnTwoNorms( Y, errorNorms );
        if( commRank == 0 )
            for( int j=0; j<numRHS; ++j )
                Output
                ("Right-hand side ",j,"\n",
                 "------------------------------------------\n",
                 "|| x     ||_2 = ",XNorms.Get(j,0),"\n",
                 "|| xComp ||_2 = ",YNorms.Get(j,0),"\n",
                 "|| A x   ||_2 = ",YOrigNorms.Get(j,0),"\n",
                 "|| error ||_2 = ",errorNorms.Get(j,0),"\n");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
