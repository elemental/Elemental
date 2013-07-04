/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Hemm.hpp"
#include "elemental/blas-like/level3/Trmm.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/blas-like/level3/TwoSidedTrmm.hpp"
#include "elemental/blas-like/level3/TwoSidedTrsm.hpp"
#include "elemental/lapack-like/HermitianGenDefiniteEig.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/matrices/Identity.hpp"
using namespace std;
using namespace elem;

template<typename F>
void TestCorrectness
( bool print,
  HermitianGenDefiniteEigType eigType,
  UpperOrLower uplo,
  const DistMatrix<F>& A,
  const DistMatrix<F>& B,
  const DistMatrix<BASE(F)>& w,
  const DistMatrix<F>& X,
  const DistMatrix<F>& AOrig,
  const DistMatrix<F>& BOrig )
{
    typedef BASE(F) R;
    const Grid& g = A.Grid();
    const int n = X.Height();
    const int k = X.Width();

    if( g.Rank() == 0 )
    {
        cout << "  Gathering computed eigenvalues...";
        cout.flush();
    }
    DistMatrix<R,MR,STAR> w_MR_STAR(true,X.RowAlignment(),g); 
    w_MR_STAR = w;
    if( g.Rank() == 0 )
        cout << "DONE" << endl;

    if( eigType == AXBX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of AX from BXW..." << endl;
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        DiagonalScale( RIGHT, NORMAL, w_MR_STAR, Y );
        // Y := Y - AX = BXW - AX
        Hemm( LEFT, uplo, F(-1), AOrig, X, F(1), Y );
        // Find the infinity norms of A, B, X, and AX-BXW
        R infNormOfA = HermitianInfinityNorm( uplo, AOrig );
        R frobNormOfA = HermitianFrobeniusNorm( uplo, AOrig );
        R infNormOfB = HermitianInfinityNorm( uplo, BOrig );
        R frobNormOfB = HermitianFrobeniusNorm( uplo, BOrig );
        R oneNormOfX = OneNorm( X );
        R infNormOfX = InfinityNorm( X );
        R frobNormOfX = FrobeniusNorm( X );
        R oneNormOfError = OneNorm( Y );
        R infNormOfError = InfinityNorm( Y );
        R frobNormOfError = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A X - B X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A X - B X W||_oo = " << infNormOfError << "\n"
                 << "    ||A X - B X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..."
                 << endl;
        }
        DistMatrix<F> Z(g);
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, F(-1), Z, F(1), Y );
        oneNormOfError = OneNorm( Y );
        infNormOfError = InfinityNorm( Y );
        frobNormOfError = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else if( eigType == ABX )
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of ABX from XW..." << endl;
        // Set Y := BX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        // Set Z := AY = ABX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), AOrig, Y, F(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const R omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const F chi = X.GetLocal(iLocal,jLocal);
                const F zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        R infNormOfA = HermitianInfinityNorm( uplo, AOrig );
        R frobNormOfA = HermitianFrobeniusNorm( uplo, AOrig );
        R infNormOfB = HermitianInfinityNorm( uplo, BOrig );
        R frobNormOfB = HermitianFrobeniusNorm( uplo, BOrig );
        R oneNormOfX = OneNorm( X );
        R infNormOfX = InfinityNorm( X );
        R frobNormOfX = FrobeniusNorm( X );
        R oneNormOfError = OneNorm( Z );
        R infNormOfError = InfinityNorm( Z );
        R frobNormOfError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||A B X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||A B X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||A B X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B..."
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, F(-1), Z, F(1), Y );
        oneNormOfError = OneNorm( Y );
        infNormOfError = InfinityNorm( Y );
        frobNormOfError = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B X - I||_F  = " << frobNormOfError << endl;
    }
    else /* eigType == BAX */
    {
        if( g.Rank() == 0 )
            cout << "  Testing for deviation of BAX from XW..." << endl;
        // Set Y := AX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
        // Set Z := BY = BAX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), BOrig, Y, F(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( int jLocal=0; jLocal<Z.LocalWidth(); ++jLocal )
        {
            const R omega = w_MR_STAR.GetLocal(jLocal,0); 
            for( int iLocal=0; iLocal<Z.LocalHeight(); ++iLocal )
            {
                const F chi = X.GetLocal(iLocal,jLocal);
                const F zeta = Z.GetLocal(iLocal,jLocal);
                Z.SetLocal(iLocal,jLocal,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        R infNormOfA = HermitianInfinityNorm( uplo, AOrig );
        R frobNormOfA = HermitianFrobeniusNorm( uplo, AOrig );
        R infNormOfB = HermitianInfinityNorm( uplo, BOrig );
        R frobNormOfB = HermitianFrobeniusNorm( uplo, BOrig );
        R oneNormOfX = OneNorm( X );
        R infNormOfX = InfinityNorm( X );
        R frobNormOfX = FrobeniusNorm( X );
        R oneNormOfError = OneNorm( Z );
        R infNormOfError = InfinityNorm( Z );
        R frobNormOfError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            cout << "    ||A||_1 = ||A||_oo = " << infNormOfA << "\n"
                 << "    ||A||_F            = " << frobNormOfA << "\n"
                 << "    ||B||_1 = ||B||_oo = " << infNormOfB << "\n"
                 << "    ||B||_F            = " << frobNormOfB << "\n"
                 << "    ||X||_1            = " << oneNormOfX << "\n"
                 << "    ||X||_oo           = " << infNormOfX << "\n"
                 << "    ||X||_F            = " << frobNormOfX << "\n"
                 << "    ||B A X - X W||_1  = " << oneNormOfError << "\n"
                 << "    ||B A X - X W||_oo = " << infNormOfError << "\n"
                 << "    ||B A X - X W||_F  = " << frobNormOfError << "\n\n"
                 << "  Testing orthonormality of eigenvectors w.r.t. B^-1..."
                 << endl;
        }
        Z = X;
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, F(-1), Z, F(1), Y );
        oneNormOfError = OneNorm( Y );
        infNormOfError = InfinityNorm( Y );
        frobNormOfError = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_1  = " << oneNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_oo = " << infNormOfError << "\n"
                 << "    ||X^H B^-1 X - I||_F  = " << frobNormOfError << endl;
    }
}

template<typename F>
void TestHermitianGenDefiniteEig
( bool testCorrectness, bool print,
  HermitianGenDefiniteEigType eigType, 
  bool onlyEigvals, UpperOrLower uplo, 
  int m, char range, BASE(F) vl, BASE(F) vu, int il, int iu, const Grid& g )
{
    typedef BASE(F) R;
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F> B(g), BOrig(g);
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> X(g);

    HermitianUniformSpectrum( A, m, 1, 10 );
    if( eigType == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( B, m, m );
        DistMatrix<F> C(g);
        Uniform( C, m, m );
        Herk( uplo, ADJOINT, F(1), C, F(0), B );
    }
    else
        HermitianUniformSpectrum( B, m, 1, 10 );

    if( testCorrectness && !onlyEigvals )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copies of original matrices...";
            cout.flush();
        }
        AOrig = A;
        BOrig = B;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting Hermitian Generalized-Definite Eigensolver...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigvals )
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, vl, vu );
    }
    else
    {
        if( range == 'A' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X );
        else if( range == 'I' )
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, il, iu );
        else
            HermitianGenDefiniteEig( eigType, uplo, A, B, w, X, vl, vu );
    }
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds." << endl;
    }
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( X, "eigenvectors:" );
    }
    if( testCorrectness && !onlyEigvals )
    {
        TestCorrectness
        ( print, eigType, uplo, A, B, w, X, AOrig, BOrig );
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--gridHeight","height of process grid",0);
        const int eigInt = Input("--eigType",
             "1 is A x = lambda B x, "
             "2 is A B x = lambda x, "
             "3 is B A x = lambda x",1);
        const bool onlyEigvals = Input 
            ("--onlyEigvals","only compute eigenvalues?",false);
        const char range = Input
            ("--range",
             "range of eigenpairs: 'A' for all, 'I' for index range, "
             "'V' for value range",'A');
        const int il = Input("--il","lower bound of index range",0);
        const int iu = Input("--iu","upper bound of index range",100);
        const double vl = Input("--vl","lower bound of value range",0.);
        const double vu = Input("--vu","upper bound of value range",100.);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const int m = Input("--height","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const Grid g( comm, r );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocal );
        SetLocalSymvBlocksize<Complex<double> >( nbLocal );
        if( range != 'A' && range != 'I' && range != 'V' )
            throw runtime_error("'range' must be 'A', 'I', or 'V'");
        if( testCorrectness && onlyEigvals && commRank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;
        HermitianGenDefiniteEigType eigType;
        string eigTypeString;
        if( eigInt == 1 )
        {
            eigType = AXBX;
            eigTypeString = "AXBX";
        }
        else if( eigInt == 2 )
        {
            eigType = ABX;
            eigTypeString = "ABX";
        }
        else if( eigInt == 3 )
        {
            eigType = BAX;
            eigTypeString = "BAX";
        }
        else
            throw logic_error("Invalid eigenvalue problem integer");
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test " 
                 << ( uplo==LOWER ? "lower" : "upper" )
                 << " " << eigTypeString << " HermitianGenDefiniteEig." << endl;

        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision normal tridiag algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestHermitianGenDefiniteEig<double>
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEig<double>
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEig<double>
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------------------\n"
                 << "Testing with double-precision complex normal algorithm:\n"
                 << "-------------------------------------------------------" 
                 << endl;
        }
        TestHermitianGenDefiniteEig<Complex<double> >
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "row-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestHermitianGenDefiniteEig<Complex<double> >
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );

        if( commRank == 0 )
        {
            cout << "---------------------------------------------------\n"
                 << "Double-precision complex square tridiag algorithm, \n"
                 << "col-major grid:\n"
                 << "---------------------------------------------------"
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestHermitianGenDefiniteEig<Complex<double> >
        ( testCorrectness, print, 
          eigType, onlyEigvals, uplo, m, range, vl, vu, il, iu, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
