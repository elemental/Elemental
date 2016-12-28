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

template<typename Field>
void TestSequentialLeastSquares
( Int numRHS, double gamma, const string& filename,
  bool feasible, bool ones, bool print )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Output("Testing with ",TypeName<Field>());

    SparseMatrix<Field> A;
    Read( A, filename, MATRIX_MARKET );
    if( print )
        Print( A, "A" );
    const Int m = A.Height();
    const Int n = A.Width();
    Output("Read matrix was ",m," x ",n," with ",A.NumEntries()," nonzeros");

    /*
    SparseMatrix<Field> ATwice;
    VCat( A, A, ATwice );
    if( print )
        Print( ATwice, "ATwice" );
    */

    Matrix<Field> B;
    if( feasible )
    {
        //Output("Generating a duplicated feasible linear system");
        Output("Generating a feasible linear system");
        Matrix<Field> X;
        if( ones )
            Ones( X, n, numRHS );
        else
            Uniform( X, n, numRHS );
        Zeros( B, m, numRHS );
        Multiply( NORMAL, Field(1), A, X, Field(0), B );
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

    /*
    Matrix<Field> BTwice;
    VCat( B, B, BTwice );
    if( print )
        Print( BTwice, "BTwice" );
    Matrix<Real> BTwiceNorms;
    ColumnTwoNorms( BTwice, BTwiceNorms );
    Print( BTwiceNorms, "BTwice column norms" );
    */
    Matrix<Real> BNorms;
    ColumnTwoNorms( B, BNorms );
    Print( BNorms, "B column norms" );

    Matrix<Field> X;
    if( gamma == double(0) )
    {
        //LeastSquares( NORMAL, ATwice, BTwice, X );
        LeastSquares( NORMAL, A, B, X );
    }
    else
    {
        //Ridge( NORMAL, ATwice, BTwice, Real(gamma), X );
        Ridge( NORMAL, A, B, Real(gamma), X );
    }
    if( print )
        Print( X, "X" );
    Matrix<Real> XNorms;
    ColumnTwoNorms( X, XNorms );
    Print( XNorms, "X norms" );

    // Compute the residual
    //Matrix<Field> E( BTwice );
    //Multiply( NORMAL, Field(-1), ATwice, X, Field(1), E );
    Matrix<Field> E( B );
    Multiply( NORMAL, Field(-1), A, X, Field(1), E );
    Matrix<Real> residNorms;
    ColumnTwoNorms( E, residNorms );
    Print( residNorms, "residual norms" );
    //DiagonalSolve( RIGHT, NORMAL, BTwiceNorms, residNorms );
    DiagonalSolve( RIGHT, NORMAL, BNorms, residNorms );
    Print( residNorms, "relative residual norms" );
    Output("Objectives:");
    for( Int j=0; j<numRHS; ++j )
      Output("  ",SafeNorm(residNorms(0,j),gamma*XNorms(0,j)));
    Output("");
}

template<typename Field>
void TestLeastSquares
( Int numRHS, double gamma, const string& filename, bool feasible, bool ones,
  bool print, const El::Grid& grid )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const int commRank = grid.Rank();
    OutputFromRoot(grid.Comm(),"Testing with ",TypeName<Field>());

    DistSparseMatrix<Field> A(grid);
    Read( A, filename, MATRIX_MARKET );
    if( print )
        Print( A, "A" );
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    OutputFromRoot
    (grid.Comm(),"Read matrix was ",m," x ",n," with ",numEntries," nonzeros");

    DistSparseMatrix<Field> ATwice(grid);
    VCat( A, A, ATwice );
    if( print )
        Print( ATwice, "ATwice" );

    DistMultiVec<Field> B(grid);
    if( feasible )
    {
        OutputFromRoot
        (grid.Comm(),"Generating a duplicated feasible linear system");
        DistMultiVec<Field> X(grid);
        if( ones )
            Ones( X, n, numRHS );
        else
            Uniform( X, n, numRHS );
        Zeros( B, m, numRHS );
        Multiply( NORMAL, Field(1), A, X, Field(0), B );
    }
    else
    {
        OutputFromRoot(grid.Comm(),"Generating a set of right-hand sides");
        if( ones )
            Ones( B, m, numRHS );
        else
            Uniform( B, m, numRHS );
    }
    if( print )
        Print( B, "B" );

    DistMultiVec<Field> BTwice(grid);
    VCat( B, B, BTwice );
    if( print )
        Print( BTwice, "BTwice" );
    Matrix<Real> BTwiceNorms;
    ColumnTwoNorms( BTwice, BTwiceNorms );
    if( commRank == 0 )
        Print( BTwiceNorms, "BTwice column norms" );

    DistMultiVec<Field> X;
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
    DistMultiVec<Field> E( BTwice );
    Multiply( NORMAL, Field(-1), ATwice, X, Field(1), E );
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
            const El::Grid grid( comm );
            TestLeastSquares<double>
            ( numRHS, gamma, filename, feasible, ones, print, grid );
#ifdef EL_HAVE_QD
            TestLeastSquares<DoubleDouble>
            ( numRHS, gamma, filename, feasible, ones, print, grid );
#ifdef EL_RELEASE
            TestLeastSquares<QuadDouble>
            ( numRHS, gamma, filename, feasible, ones, print, grid );
#endif
#endif
#ifdef EL_HAVE_QUAD
#ifdef EL_RELEASE
            TestLeastSquares<Quad>
            ( numRHS, gamma, filename, feasible, ones, print, grid );
#endif
#endif
#ifdef EL_HAVE_MPC
#ifdef EL_RELEASE
            TestLeastSquares<BigFloat>
            ( numRHS, gamma, filename, feasible, ones, print, grid );
#endif
#endif
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
