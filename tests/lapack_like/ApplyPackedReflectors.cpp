/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void TestCorrectness
( LeftOrRight side,
  UpperOrLower uplo,
  ForwardOrBackward order,
  Conjugation conjugation,
  Int offset,
  bool printMatrices,
  const DistMatrix<F>& H,
  const DistMatrix<F,MD,STAR>& householderScalars )
{
    typedef Base<F> Real;
    const Grid& g = H.Grid();
    const Int m = H.Height();
    const Real eps = limits::Epsilon<Real>();

    OutputFromRoot(g.Comm(),"Testing orthogonality of transform...");
    PushIndent();

    // Form Z := Q^H Q or Q Q^H as an approximation to identity
    DistMatrix<F> Y(g);
    Identity( Y, m, m );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, order, conjugation, offset,
      H, householderScalars, Y );
    if( printMatrices )
    {
        DistMatrix<F> W(g);
        Identity( W, m, m );
        if( order == FORWARD )
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, BACKWARD, conjugation, offset,
              H, householderScalars, W );
            Print( Y, "Q" );
            Print( W, "Q^H" );
        }
        else
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, FORWARD, conjugation, offset,
              H, householderScalars, W );
            Print( Y, "Q^H" );
            Print( W, "Q" );
        }
    }
    DistMatrix<F> Z(g);
    Zeros( Z, m, m );
    Herk( uplo, NORMAL, Real(1), Y, Real(0), Z );
    MakeHermitian( uplo, Z );

    // Form X := -I + Q^H Q or Q Q^H
    ShiftDiagonal( Z, F(-1) );
    if( printMatrices )
    {
        if( order == FORWARD )
            Print( Z, "Q Q^H - I" );
        else
            Print( Z, "Q^H Q - I" );
    }

    // Compute the maximum deviance
    const Real infNormError = InfinityNorm( Z );
    const Real relError = infNormError / (m*eps);
    if( order == FORWARD )
    {
        OutputFromRoot(g.Comm(), "||Q Q^H - I||_oo / (m eps) = ",relError);
    }
    else
    {
        OutputFromRoot(g.Comm(), "||Q^H Q - I||_oo / (m eps) = ",relError);
    }
    PopIndent();
    // TODO: Use a more refined failure condition
    if( relError > Real(10) )
        LogicError("Relative error was unacceptably large");
}

template<typename F>
void TestUT
( const Grid& g,
  LeftOrRight side,
  UpperOrLower uplo,
  ForwardOrBackward order,
  Conjugation conjugation,
  Int m,
  Int offset,
  bool correctness,
  bool printMatrices )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();
    DistMatrix<F> H(g), A(g);
    Uniform( H, m, m );
    Uniform( A, m, m );

    const Int diagLength = DiagonalLength(H.Height(),H.Width(),offset);
    DistMatrix<F,MD,STAR> householderScalars(g);
    householderScalars.SetRoot( H.DiagonalRoot(offset) );
    householderScalars.AlignCols( H.DiagonalAlign(offset) );
    householderScalars.Resize( diagLength, 1 );

    DistMatrix<F> HCol(g);
    if( uplo == LOWER )
    {
        for( Int i=0; i<householderScalars.Height(); ++i )
        {
            // View below the diagonal containing the implicit 1
            HCol = View( H, i-offset+1, i, m-(i-offset+1), 1 );
            F norm = Nrm2( HCol );
            F alpha = F(2)/(norm*norm+F(1));
            householderScalars.Set( i, 0, alpha );
        }
    }
    else
    {
        for( Int i=0; i<householderScalars.Height(); ++i )
        {
            // View above the diagonal containing the implicit 1
            HCol = View( H, 0, i+offset, i, 1 );
            F norm = Nrm2( HCol );
            F alpha = F(2)/(norm*norm+F(1));
            householderScalars.Set( i, 0, alpha );
        }
    }

    if( printMatrices )
    {
        Print( H, "H" );
        Print( A, "A" );
        Print( householderScalars, "householderScalars" );
    }

    OutputFromRoot(g.Comm(),"Starting UT transform...");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, order, conjugation, offset,
      H, householderScalars, A );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = 8.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds (",gFlops," GFlop/s)");
    if( printMatrices )
        Print( A, "A after factorization" );
    if( correctness )
    {
        TestCorrectness
        ( side, uplo, order, conjugation, offset, printMatrices, H, householderScalars );
    }
    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        Int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","store in triangle: L/U",'L');
        const bool forward = Input("--forward","forward application?",true);
        const bool conjugate = Input("--conjugate","conjugate?",false);
        const Int m = Input("--height","height of matrix",100);
        const Int offset = Input("--offset","diagonal offset for storage",0);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool correctness  = Input
            ("--correctness","test correctness?",true);
        const bool printMatrices = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
#endif

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const ForwardOrBackward dir = ( forward ? FORWARD : BACKWARD );
        const Conjugation conjugation =
            ( conjugate ? CONJUGATED : UNCONJUGATED );
        SetBlocksize( nb );
        if( uplo == LOWER && offset > 0 )
            LogicError
            ("Offset cannot be positive if transforms are in lower triangle");
        else if( uplo == UPPER && offset < 0 )
            LogicError
            ("Offset cannot be negative if transforms are in upper triangle");

        ComplainIfDebug();
        OutputFromRoot(g.Comm(),"Will test UT transform");

        TestUT<float>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<Complex<float>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );

        TestUT<double>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<Complex<double>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );

#ifdef EL_HAVE_QD
        TestUT<DoubleDouble>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<QuadDouble>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );

        TestUT<Complex<DoubleDouble>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<Complex<QuadDouble>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
#endif

#ifdef EL_HAVE_QUAD
        TestUT<Quad>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<Complex<Quad>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
#endif

#ifdef EL_HAVE_MPC
        TestUT<BigFloat>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
        TestUT<Complex<BigFloat>>
        ( g, side, uplo, dir, conjugation, m, offset,
          correctness, printMatrices );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
