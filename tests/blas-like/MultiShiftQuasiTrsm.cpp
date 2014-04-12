/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_MAKETRAPEZOIDAL_INC
#include ELEM_GEMM_INC
#include ELEM_MULTISHIFTQUASITRSM_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
#include ELEM_UNIFORM_INC
using namespace std;
using namespace elem;

template<typename F> 
void TestMultiShiftQuasiTrsm
( bool print,
  LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  Int m, Int n, F alpha, const Grid& g )
{
    typedef Base<F> Real;
    DistMatrix<F> A(g), X(g);
    DistMatrix<F,VR,STAR> shifts(g);

    if( side == LEFT )
    {
        HermitianUniformSpectrum( A, m, 1, 10 );
        Uniform( shifts, m, 1, F(0), Real(0.5) );
    }
    else
    {
        HermitianUniformSpectrum( A, n, 1, 10 );
        Uniform( shifts, n, 1, F(0), Real(0.5) );
    }
    // TODO: Enforce fact that A should not have consecutive nonzeros in the
    //       subdiagonal
    auto H( A );
    if( uplo == LOWER )
        MakeTrapezoidal( LOWER, H, 1 );
    else
        MakeTrapezoidal( UPPER, H, -1 );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    Zeros( Y, m, n );
    if( side == LEFT )
    {
        Gemm( NORMAL, NORMAL, F(1), H, X, F(1), Y );
        for( Int j=0; j<n; ++j )
        {
            auto x = LockedView( X, 0, j, m, 1 );
            auto y =       View( Y, 0, j, m, 1 );
            Axpy( -shifts.Get(j,0), x, y );
        }
    }
    else
    {
        Gemm( NORMAL, NORMAL, F(1), X, H, F(1), Y );
        for( Int i=0; i<m; ++i )
        {
            auto x = LockedView( X, i, 0, 1, n );
            auto y =       View( Y, i, 0, 1, n );
            Axpy( -shifts.Get(i,0), x, y );
        }
    }

    if( print )
    {
        Print( A, "A" );
        Print( H, "H" );
        Print( shifts, "shifts" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting MultiShiftQuasiTrsm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    MultiShiftQuasiTrsm( side, uplo, orientation, alpha, A, shifts, Y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 
        ( side==LEFT ? double(m)*double(m)*double(n)
                     : double(m)*double(n)*double(n) ) /(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. \n"
             << "  Time = " << runTime << " seconds. GFlops ~= " << gFlops 
             << endl;
    }
    if( print )
        Print( Y, "Y after solve" );
    Axpy( F(-1), X, Y );
    const auto HFrob = FrobeniusNorm( H );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        cout << "|| H ||_F = " << HFrob << "\n"
             << "|| X ||_F = " << XFrob << "\n"
             << "|| E ||_F = " << EFrob << "\n" << std::endl;
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char sideChar = Input("--side","side to solve from: L/R",'L');
        const char uploChar = Input
            ("--uplo","lower or upper quasi-triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of quasi-triangular matrix: N/T/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test MultiShiftQuasiTrsm" 
                 << sideChar << uploChar << transChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestMultiShiftQuasiTrsm<double>
        ( print, side, uplo, orientation, m, n, 3., g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestMultiShiftQuasiTrsm<Complex<double>>
        ( print, side, uplo, orientation, m, n, Complex<double>(3), g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
