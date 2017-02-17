/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Field>
void TestCorrectness
( UpperOrLower uplo,
  UnitOrNonUnit diag,
  const Matrix<Field>& A,
  const Matrix<Field>& AOrig,
  bool print )
{
    typedef Base<Field> Real;
    const Int m = AOrig.Height();
    const Real oneNormA = OneNorm( AOrig );
    const Real eps = limits::Epsilon<Real>();

    // Test I - inv(A) A
    Matrix<Field> X;
    X = AOrig;
    MakeTrapezoidal( uplo, X );
    Trmm( LEFT, uplo, NORMAL, diag, Field(1), A, X );
    ShiftDiagonal( X, Field(-1) );

    const Real maxError = MaxNorm( X );
    const Real relError = maxError / (eps*m*oneNormA);
    Output("||I - inv(A) A||_max / (eps m ||A||_1) = ",relError);

    // TODO: More rigorous failure condition
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestCorrectness
( UpperOrLower uplo,
  UnitOrNonUnit diag,
  const DistMatrix<Field>& A,
  const DistMatrix<Field>& AOrig,
  bool print )
{
    typedef Base<Field> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();
    const Real oneNormA = OneNorm( AOrig );
    const Real eps = limits::Epsilon<Real>();

    // Test I - inv(A) A
    DistMatrix<Field> X(g);
    X = AOrig;
    MakeTrapezoidal( uplo, X );
    Trmm( LEFT, uplo, NORMAL, diag, Field(1), A, X );
    ShiftDiagonal( X, Field(-1) );

    const Real maxError = MaxNorm( X );
    const Real relError = maxError / (eps*m*oneNormA);
    OutputFromRoot
    (g.Comm(),
     "||I - inv(A) A||_max / (eps m ||A||_1) = ",relError);

    // TODO: More rigorous failure condition
    if( relError > Real(10) )
        LogicError("Unacceptably large relative error");
}

template<typename Field>
void TestTriangularInverse
( UpperOrLower uplo,
  UnitOrNonUnit diag,
  Int m,
  bool correctness,
  bool print )
{
    Output("Testing with ",TypeName<Field>());
    PushIndent();

    Matrix<Field> A, AOrig;
    Uniform( A, m, m );
    MakeTrapezoidal( uplo, A );
    ShiftDiagonal( A, Field(3) );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    Output("Starting triangular inversion...");
    Timer timer;
    timer.Start();
    TriangularInverse( uplo, diag, A );
    const double runTime = timer.Stop();
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    Output("Time = ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( A, "A after inversion" );
    if( correctness )
        TestCorrectness( uplo, diag, A, AOrig, print );
    PopIndent();
}

template<typename Field>
void TestTriangularInverse
( const Grid& g,
  UpperOrLower uplo,
  UnitOrNonUnit diag,
  Int m,
  bool correctness,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<Field>());
    PushIndent();

    DistMatrix<Field> A(g), AOrig(g);
    Uniform( A, m, m );
    MakeTrapezoidal( uplo, A );
    ShiftDiagonal( A, Field(3) );
    if( correctness )
        AOrig = A;
    if( print )
        Print( A, "A" );

    OutputFromRoot(g.Comm(),"Starting triangular inversion...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    TriangularInverse( uplo, diag, A );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = IsComplex<Field>::value ? 4*realGFlops : realGFlops;
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( A, "A after inversion" );
    if( correctness )
        TestCorrectness( uplo, diag, A, AOrig, print );
    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool correctness =
          Input("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        OutputFromRoot
        (g.Comm(),"Will test TriangularInverse",uploChar,diagChar);

        if( sequential && mpi::Rank() == 0 )
        {
            TestTriangularInverse<float>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<Complex<float>>
            ( uplo, diag, m, correctness, print );

            TestTriangularInverse<double>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<Complex<double>>
            ( uplo, diag, m, correctness, print );

#ifdef EL_HAVE_QD
            TestTriangularInverse<DoubleDouble>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<QuadDouble>
            ( uplo, diag, m, correctness, print );

            TestTriangularInverse<Complex<DoubleDouble>>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<Complex<QuadDouble>>
            ( uplo, diag, m, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
            TestTriangularInverse<Quad>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<Complex<Quad>>
            ( uplo, diag, m, correctness, print );
#endif

#ifdef EL_HAVE_MPC
            TestTriangularInverse<BigFloat>
            ( uplo, diag, m, correctness, print );
            TestTriangularInverse<Complex<BigFloat>>
            ( uplo, diag, m, correctness, print );
#endif
        }

        TestTriangularInverse<float>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<Complex<float>>
        ( g, uplo, diag, m, correctness, print );

        TestTriangularInverse<double>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<Complex<double>>
        ( g, uplo, diag, m, correctness, print );

#ifdef EL_HAVE_QD
        TestTriangularInverse<DoubleDouble>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<QuadDouble>
        ( g, uplo, diag, m, correctness, print );

        TestTriangularInverse<Complex<DoubleDouble>>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<Complex<QuadDouble>>
        ( g, uplo, diag, m, correctness, print );
#endif

#ifdef EL_HAVE_QUAD
        TestTriangularInverse<Quad>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<Complex<Quad>>
        ( g, uplo, diag, m, correctness, print );
#endif

#ifdef EL_HAVE_MPC
        TestTriangularInverse<BigFloat>
        ( g, uplo, diag, m, correctness, print );
        TestTriangularInverse<Complex<BigFloat>>
        ( g, uplo, diag, m, correctness, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
