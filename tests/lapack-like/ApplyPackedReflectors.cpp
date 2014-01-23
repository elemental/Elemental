/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level1/Nrm2.hpp"
#include "elemental/blas-like/level1/UpdateDiagonal.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( LeftOrRight side, UpperOrLower uplo, ForwardOrBackward order,
  Conjugation conjugation, Int offset, bool printMatrices,
  const DistMatrix<F>& H,
  const DistMatrix<F,MD,STAR>& t )
{
    typedef Base<F> Real;
    const Grid& g = H.Grid();
    const Int m = H.Height();

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of transform..." << endl;

    // Form Z := Q^H Q or Q Q^H as an approximation to identity
    auto Y = Identity<F>( g, m, m );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, order, conjugation, offset, H, t, Y );
    if( printMatrices )
    {
        auto W = Identity<F>( g, m, m );
        if( order == FORWARD )
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, BACKWARD, conjugation, offset, H, t, W );
            Print( Y, "Q" );
            Print( W, "Q^H" );
        }
        else
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, FORWARD, conjugation, offset, H, t, W );
            Print( Y, "Q^H" );
            Print( W, "Q" );
        }
    }
    auto Z = Zeros<F>( g, m, m );
    Herk( uplo, NORMAL, F(1), Y, F(0), Z );
    MakeHermitian( uplo, Z );
    
    // Form X := -I + Q^H Q or Q Q^H
    UpdateDiagonal( Z, F(-1) );
    if( printMatrices )
    {
        if( order == FORWARD )
            Print( Z, "Q Q^H - I" );
        else
            Print( Z, "Q^H Q - I" );
    }

    // Compute the maximum deviance
    const Real oneNormOfError = OneNorm( Z );
    const Real infNormOfError = InfinityNorm( Z );
    const Real frobNormOfError = FrobeniusNorm( Z );
    if( g.Rank() == 0 )
    {
        if( order == FORWARD )
        {
            cout << "    ||Q Q^H - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q Q^H - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q Q^H - I||_F  = " << frobNormOfError << endl;
        }
        else
        {
            cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
                 << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
                 << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
        }
    }
}

template<typename F>
void TestUT
( LeftOrRight side, UpperOrLower uplo, 
  ForwardOrBackward order, Conjugation conjugation,
  Int m, Int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    auto H = Uniform<F>( g, m, m );
    auto A = Uniform<F>( g, m, m );

    const Int diagLength = DiagonalLength(H.Height(),H.Width(),offset);
    DistMatrix<F,MD,STAR> t(g);
    t.AlignWithDiagonal( H, offset );
    t.Resize( diagLength, 1 );

    DistMatrix<F> HCol(g);
    if( uplo == LOWER )
    {
        for( Int i=0; i<t.Height(); ++i )
        {
            // View below the diagonal containing the implicit 1
            HCol = View( H, i-offset+1, i, m-(i-offset+1), 1 );
            F norm = Nrm2( HCol );
            F alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }
    else
    {
        for( Int i=0; i<t.Height(); ++i ) 
        {
            // View above the diagonal containing the implicit 1
            HCol = View( H, 0, i+offset, i, 1 );
            F norm = Nrm2( HCol );
            F alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }

    if( printMatrices )
    {
        Print( H, "H" );
        Print( A, "A" );
        Print( t, "t" );
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, order, conjugation, offset, H, t, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 8.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        Print( A, "A after factorization" );
    if( testCorrectness )
    {
        TestCorrectness
        ( side, uplo, order, conjugation, offset, printMatrices, H, t );
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );
    const Int commSize = mpi::CommSize( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","store in triangle: L/U",'L');
        const bool forward = Input("--forward","forward application?",true);
        const bool conjugate = Input("--conjugate","conjugate?",false);
        const Int m = Input("--height","height of matrix",100);
        const Int offset = Input("--offset","diagonal offset for storage",0);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness  = Input
            ("--correctness","test correctness?",true);
        const bool printMatrices = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const Grid g( comm, r );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const ForwardOrBackward order = ( forward ? FORWARD : BACKWARD );
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
        if( commRank == 0 )
            cout << "Will test UT transform" << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestUT<double>
        ( side, uplo, order, conjugation, m, offset, 
          testCorrectness, printMatrices, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestUT<Complex<double>>
        ( side, uplo, order, conjugation, m, offset, 
          testCorrectness, printMatrices, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
