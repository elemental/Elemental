/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Test if (op(H) - mu_j I) x_j = y_j for each j.
// This is checked by testing the norm of  op(H) X - X Mu - Y.
template<typename F> 
void TestCorrectness
( UpperOrLower uplo,
  Orientation orientation,
  const DistMatrix<F,VC,STAR>& H, 
  const DistMatrix<F,VR,STAR>& shifts, 
  const DistMatrix<F,STAR,VR>& X, 
  const DistMatrix<F,STAR,VR>& Y,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();

    auto modShifts( shifts );
    if( orientation == ADJOINT )
        Conjugate( modShifts );
    
    DistMatrix<F> Z( Y );
    for( Int j=0; j<n; ++j )
    {
        auto x = LockedView( X, 0, j, m, 1 );
        auto z =       View( Z, 0, j, m, 1 );
        Axpy( modShifts.Get(j,0), x, z );
    }
    {
        DistMatrix<F> H_MC_MR( H ), X_MC_MR(X);
        Gemm( orientation, NORMAL, F(-1), H_MC_MR, X_MC_MR, F(1), Z );
    }

    if( print )
    {
        Print( H, "H" );
        Print( X, "X" );
        Print( Y, "Y" );
        Print( shifts, "shifts" );
        Print( Z, "-H X + X Mu + Y" );
    }
    if( display )
    {
        Display( H, "H" );
        Display( X, "X" );
        Display( Y, "Y" );
        Display( shifts, "shifts" );
        Display( Z, "-H X + X Mu + Y" );
    }

    const Real YFrob = FrobeniusNorm( Y );
    const Real YInf = InfinityNorm( Y );
    const Real HFrob = FrobeniusNorm( H );
    const Real HInf = InfinityNorm( H );
    const Real ZFrob = FrobeniusNorm( Z );
    const Real ZInf = InfinityNorm( Z );
    if( mpi::Rank() == 0 )
        Output
        ("    || H ||_F  = ",HFrob,"\n",
         "    || H ||_oo = ",HInf,"\n",
         "    || Y ||_F  = ",YFrob,"\n",
         "    || Y ||_oo = ",YInf,"\n",
         "    || H X - X Mu - Y ||_F  = ",ZFrob,"\n",
         "    || H X - X Mu - Y ||_oo = ",ZInf);
}

template<typename F>
void TestHessenberg
( const Grid& g,
  UpperOrLower uplo,
  Orientation orientation,
  Int m,
  Int n, 
  bool testCorrectness,
  bool print,
  bool display )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    DistMatrix<F,VC,STAR> H(g);
    DistMatrix<F,STAR,VR> X(g), Y(g);
    DistMatrix<F,VR,STAR> shifts(g);

    Uniform( H, m, m );
    ShiftDiagonal( H, F(5) ); // ensure that H-mu is far from zero
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );

    Uniform( X, m, n );
    Uniform( Y, m, n );
    Uniform( shifts, n, 1 );

    X = Y;
    if( mpi::Rank() == 0 )
        Output("  Starting Hessenberg solve...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    MultiShiftHessSolve( uplo, orientation, F(1), H, shifts, X );
    mpi::Barrier( mpi::COMM_WORLD );
    const double runTime = mpi::Time() - startTime;
    // TODO: Flop calculation
    if( mpi::Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( testCorrectness )
        TestCorrectness( uplo, orientation, H, shifts, X, Y, print, display );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char orientChar = Input("--orient","orientation: N/T/C",'N');
        const Int m = Input("--m","height of Hessenberg matrix",100);
        const Int n = Input("--n","number of right-hand sides",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid grid( comm, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orient = CharToOrientation( orientChar );
        SetBlocksize( nb );
        ComplainIfDebug();

        TestHessenberg<float>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
        TestHessenberg<Complex<float>>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );

        TestHessenberg<double>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
        TestHessenberg<Complex<double>>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );

#ifdef EL_HAVE_QD
        TestHessenberg<DoubleDouble>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
        TestHessenberg<QuadDouble>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_QUAD
        TestHessenberg<Quad>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
        TestHessenberg<Complex<Quad>>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_MPC
        TestHessenberg<BigFloat>
        ( grid, uplo, orient, m, n, testCorrectness, print, display );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
