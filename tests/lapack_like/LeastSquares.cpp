/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

// TODO(poulson): Add ability to tune reg0Tmp, reg0Perm, reg1Tmp, etc. within
// LeastSquaresCtrl

template<typename F>
void TestSequentialLeastSquares
( Int numRHS, double gamma, const string& filename,
  bool feasible, bool ones, bool print )
{
    DEBUG_CSE
    typedef Base<F> Real;
    Output("Testing with ",TypeName<F>());

    SparseMatrix<F> A;
    Read( A, filename, MATRIX_MARKET );
    if( print )
        Print( A, "A" );
    const Int m = A.Height();
    const Int n = A.Width();
    Output("Read matrix was ",m," x ",n," with ",A.NumEntries()," nonzeros");

    SparseMatrix<F> ATwice;
    VCat( A, A, ATwice );
    if( print )
        Print( ATwice, "ATwice" );

    Matrix<F> B;
    if( feasible )
    {
        Output("Generating a duplicated feasible linear system");
        Matrix<F> X;
        if( ones )
            Ones( X, n, numRHS );
        else
            Uniform( X, n, numRHS );
        Zeros( B, m, numRHS );
        Multiply( NORMAL, F(1), A, X, F(0), B );
    }
    else
    {
        Output("Generating a set of right-hand sides");
        if( ones )
            Ones( B, m, numRHS );
        else
            Uniform( B, m, numRHS );
    }
    if( print )
        Print( B, "B" );

    Matrix<F> BTwice;
    VCat( B, B, BTwice );
    if( print )
        Print( BTwice, "BTwice" );
    Matrix<Real> BTwiceNorms;
    ColumnTwoNorms( BTwice, BTwiceNorms );
    Print( BTwiceNorms, "BTwice column norms" );

    Matrix<F> X;
    if( gamma == double(0) )
    {
        LeastSquares( NORMAL, ATwice, BTwice, X );
    }
    else
    {
        Ridge( NORMAL, ATwice, BTwice, Real(gamma), X );
    }
    if( print )
        Print( X, "X" );
    Matrix<Real> XNorms;
    ColumnTwoNorms( X, XNorms );
    Print( XNorms, "X norms" );

    // Compute the residual
    Matrix<F> E( BTwice );
    Multiply( NORMAL, F(-1), ATwice, X, F(1), E ); 
    Matrix<Real> residNorms;
    ColumnTwoNorms( E, residNorms );
    Print( residNorms, "residual norms" );
    DiagonalSolve( RIGHT, NORMAL, BTwiceNorms, residNorms );
    Print( residNorms, "relative residual norms" );
    Output("Objectives:");
    for( Int j=0; j<numRHS; ++j )
      Output("  ",SafeNorm(residNorms(0,j),gamma*XNorms(0,j)));
    Output("");
}

template<typename F>
void TestLeastSquares
( Int numRHS, double gamma, const string& filename, bool feasible, bool ones,
  bool print, mpi::Comm comm )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const int commRank = mpi::Rank(comm);
    OutputFromRoot(comm,"Testing with ",TypeName<F>());

    DistSparseMatrix<F> A(comm);
    Read( A, filename, MATRIX_MARKET );
    if( print )
        Print( A, "A" );
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    OutputFromRoot
    (comm,"Read matrix was ",m," x ",n," with ",numEntries," nonzeros");

    DistSparseMatrix<F> ATwice(comm);
    VCat( A, A, ATwice );
    if( print )
        Print( ATwice, "ATwice" );

    DistMultiVec<F> B(comm);
    if( feasible )
    {
        OutputFromRoot(comm,"Generating a duplicated feasible linear system");
        DistMultiVec<F> X(comm);
        if( ones )
            Ones( X, n, numRHS );
        else
            Uniform( X, n, numRHS );
        Zeros( B, m, numRHS );
        Multiply( NORMAL, F(1), A, X, F(0), B );
    }
    else
    {
        OutputFromRoot(comm,"Generating a set of right-hand sides");
        if( ones )
            Ones( B, m, numRHS );
        else
            Uniform( B, m, numRHS );
    }
    if( print )
        Print( B, "B" );

    DistMultiVec<F> BTwice(comm);
    VCat( B, B, BTwice );
    if( print )
        Print( BTwice, "BTwice" );
    Matrix<Real> BTwiceNorms;
    ColumnTwoNorms( BTwice, BTwiceNorms );
    if( commRank == 0 )
        Print( BTwiceNorms, "BTwice column norms" );

    DistMultiVec<F> X;
    if( gamma == double(0) )
    {
        LeastSquares( NORMAL, ATwice, BTwice, X );
    }
    else
    {
        Ridge( NORMAL, ATwice, BTwice, Real(gamma), X );
    }
    if( print )
        Print( X, "X" );
    Matrix<Real> XNorms;
    ColumnTwoNorms( X, XNorms );
    if( commRank == 0 )
        Print( XNorms, "X norms" );

    // Compute the objectives,
    //
    //   || [A; gamma*I] x - [b; 0] ||_2 =
    //     sqrt( || A x - b ||_2^2 + gamma^2 || x ||_2^2 ).
    //
    DistMultiVec<F> E( BTwice );
    Multiply( NORMAL, F(-1), ATwice, X, F(1), E ); 
    Matrix<Real> residNorms;
    ColumnTwoNorms( E, residNorms );
    if( commRank == 0 )
    {
        Print( residNorms, "residual norms" );
        DiagonalSolve( RIGHT, NORMAL, BTwiceNorms, residNorms );
        Print( residNorms, "relative residual norms" );

        Output("Objectives:");
        for( Int j=0; j<numRHS; ++j )
          Output("  ",SafeNorm(residNorms(0,j),gamma*XNorms(0,j)));
        Output("");
    }
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank(comm);

    try
    {
        const Int numRHS = Input("--numRHS","num RHS",1);
        const double gamma = Input("--gamma","regularization",0.001); 
        const string filename =
          Input
          ("--filename","path to Matrix Market",
           string("../data/lapack_like/c-41.mtx"));
        const bool feasible =
          Input("--feasible","generate a feasible RHS?",true);
        const bool ones =
          Input("--ones","use all ones for generating RHS?",false);
        const bool print = Input("--print","print matrices?",false);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool distributed =
          Input("--distributed","test distributed?",true);
        ProcessInput();

        if( sequential && commRank == 0 )
        {
            TestSequentialLeastSquares<double>
            ( numRHS, gamma, filename, feasible, ones, print );
#ifdef EL_HAVE_QD
            // We will be humane and only run up to DoubleDouble in Debug builds due to
            // excessive runtimes
            TestSequentialLeastSquares<DoubleDouble>
            ( numRHS, gamma, filename, feasible, ones, print );
#ifdef EL_RELEASE 
            TestSequentialLeastSquares<QuadDouble>
            ( numRHS, gamma, filename, feasible, ones, print );
#endif
#endif
#ifdef EL_HAVE_QUAD
#ifdef EL_RELEASE
            TestSequentialLeastSquares<Quad>
            ( numRHS, gamma, filename, feasible, ones, print );
#endif
#endif
#ifdef EL_HAVE_MPC
#ifdef EL_RELEASE
            TestSequentialLeastSquares<BigFloat>
            ( numRHS, gamma, filename, feasible, ones, print );
#endif
#endif
        }
        if( distributed )
        {
            TestLeastSquares<double>
            ( numRHS, gamma, filename, feasible, ones, print, comm );
#ifdef EL_HAVE_QD
            TestLeastSquares<DoubleDouble>
            ( numRHS, gamma, filename, feasible, ones, print, comm );
#ifdef EL_RELEASE
            TestLeastSquares<QuadDouble>
            ( numRHS, gamma, filename, feasible, ones, print, comm );
#endif
#endif
#ifdef EL_HAVE_QUAD
#ifdef EL_RELEASE
            TestLeastSquares<Quad>
            ( numRHS, gamma, filename, feasible, ones, print, comm );
#endif
#endif
#ifdef EL_HAVE_MPC
#ifdef EL_RELEASE
            TestLeastSquares<BigFloat>
            ( numRHS, gamma, filename, feasible, ones, print, comm );
#endif
#endif
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
