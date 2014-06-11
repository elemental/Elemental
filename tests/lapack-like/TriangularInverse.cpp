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
void TestCorrectness
( bool print,
  UpperOrLower uplo, UnitOrNonUnit diag,
  const DistMatrix<F>& A,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    DistMatrix<F> X(g), Y(g);
    Uniform( X, m, 100 );
    Y = X;

    // Since A o A^-1 = I, test the change introduced by the approximate comp.
    Trmm( LEFT, uplo, NORMAL, diag, F(1), A,     Y );
    Trmm( LEFT, uplo, NORMAL, diag, F(1), AOrig, Y );
    Axpy( F(-1), X, Y );

    const Real oneNormOrig = OneNorm( AOrig );
    const Real infNormOrig = InfinityNorm( AOrig );
    const Real frobNormOrig = FrobeniusNorm( AOrig );
    const Real oneNormFinal = OneNorm( A );
    const Real infNormFinal = InfinityNorm( A );
    const Real frobNormFinal = FrobeniusNorm( A );
    const Real oneNormOfError = OneNorm( Y );
    const Real infNormOfError = InfinityNorm( Y );
    const Real frobNormOfError = FrobeniusNorm( Y );
    if( g.Rank() == 0 )
    {
        cout << "||A||_1           = " << oneNormOrig << "\n"
             << "||A||_oo          = " << infNormOrig << "\n"
             << "||A||_F           = " << frobNormOrig << "\n"
             << "||A^-1||_1        = " << oneNormFinal << "\n"
             << "||A^-1||_oo       = " << infNormFinal << "\n"
             << "||A^-1||_F        = " << frobNormFinal << "\n"
             << "||A A^-1 - I||_1  = " << oneNormOfError << "\n"
             << "||A A^-1 - I||_oo = " << infNormOfError << "\n"
             << "||A A^-1 - I||_F  = " << frobNormOfError << endl;
    }
}

template<typename F> 
void TestTriangularInverse
( bool testCorrectness, bool print,
  UpperOrLower uplo, UnitOrNonUnit diag, Int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    HermitianUniformSpectrum( A, m, 1, 10 );
    MakeTriangular( uplo, A );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
    {
        cout << "  Starting triangular inversion...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TriangularInverse( uplo, diag, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 1./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        Print( A, "A after inversion" );
    if( testCorrectness )
        TestCorrectness( print, uplo, diag, A, AOrig );
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
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test TriangularInverse" << uploChar << diagChar 
                 << endl;

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestTriangularInverse<double>
        ( testCorrectness, print, uplo, diag, m, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestTriangularInverse<Complex<double>>
        ( testCorrectness, print, uplo, diag, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
