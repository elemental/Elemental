/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_HERMITIANUNIFORMSPECTRUM_INC

using namespace std;
using namespace El;

template<typename F> 
void TestMultiShiftTrsm
( bool print,
  LeftOrRight side, UpperOrLower uplo, Orientation orientation, 
  Int m, Int n, F alpha, const Grid& g )
{
    typedef Base<F> Real;
    DistMatrix<F> U(g), X(g);
    DistMatrix<F,VR,STAR> shifts(g);

    if( side == LEFT )
    {
        HermitianUniformSpectrum( U, m, 1, 10 );
        Uniform( shifts, n, 1, F(0), Real(0.5) );
    }
    else
    {
        HermitianUniformSpectrum( U, n, 1, 10 );
        Uniform( shifts, m, 1, F(0), Real(0.5) );
    }
    MakeTriangular( uplo, U );

    auto modShifts(shifts);
    if( orientation == ADJOINT )
        Conjugate( modShifts );

    Uniform( X, m, n );
    DistMatrix<F> Y(g);
    Zeros( Y, m, n );
    if( side == LEFT )
    {
        Gemm( orientation, NORMAL, F(1)/alpha, U, X, F(1), Y );
        for( Int j=0; j<n; ++j )
        {
            auto x = LockedView( X, 0, j, m, 1 );
            auto y =       View( Y, 0, j, m, 1 );
            Axpy( -modShifts.Get(j,0)/alpha, x, y );
        }
    }
    else
    {
        Gemm( NORMAL, orientation, F(1)/alpha, X, U, F(1), Y );
        for( Int i=0; i<m; ++i )
        {
            auto x = LockedView( X, i, 0, 1, n );
            auto y =       View( Y, i, 0, 1, n );
            Axpy( -modShifts.Get(i,0)/alpha, x, y );
        }
    }

    if( print )
    {
        Print( U, "U" );
        Print( shifts, "shifts" );
        Print( X, "X" );
        Print( Y, "Y" );
    }
    if( g.Rank() == 0 )
    {
        cout << "  Starting MultiShiftTrsm...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    MultiShiftTrsm( side, uplo, orientation, alpha, U, shifts, Y );
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
    const auto UFrob = FrobeniusNorm( U );
    const auto XFrob = FrobeniusNorm( X );
    const auto EFrob = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        cout << "|| U ||_F = " << UFrob << "\n"
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
            cout << "Will test MultiShiftTrsm" 
                 << sideChar << uploChar << transChar << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestMultiShiftTrsm<double>( print, side, uplo, orientation, m, n, 3., g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestMultiShiftTrsm<Complex<double>>
        ( print, side, uplo, orientation, m, n, Complex<double>(3), g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
