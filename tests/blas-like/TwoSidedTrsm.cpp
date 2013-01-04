/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( bool print, UpperOrLower uplo, UnitOrNonUnit diag,
  const DistMatrix<F>& A, const DistMatrix<F>& B, const DistMatrix<F>& AOrig )
{
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    const int k=100;
    DistMatrix<F> X(g), Y(g), Z(g);
    Uniform( m, k, X );
    Y = X;
    Zeros( m, k, Z );

    if( uplo == LOWER )
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // tril(B)^-1 AOrig tril(B)^-H
        Trsm( LEFT, LOWER, ADJOINT, diag, F(1), B, Y );
        Hemm( LEFT, LOWER, F(1), AOrig, Y, F(0), Z );
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), B, Z );
        Hemm( LEFT, LOWER, F(-1), A, X, F(1), Z );
        R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        R infNormOfA = HermitianNorm( uplo, A, INFINITY_NORM );
        R frobNormOfA = HermitianNorm( uplo, A, FROBENIUS_NORM );
        R oneNormOfError = Norm( Z, ONE_NORM );
        R infNormOfError = Norm( Z, INFINITY_NORM );
        R frobNormOfError = Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "||AOrig||_1 = ||AOrig||_oo     = "
                 << infNormOfAOrig << "\n"
                 << "||AOrig||_F                    = "
                 << frobNormOfAOrig << "\n"
                 << "||A||_1 = ||A||_oo             = "
                 << infNormOfA << "\n"
                 << "||A||_F                        = "
                 << frobNormOfA << "\n"
                 << "||A X - L^-1 AOrig L^-H X||_1  = "
                 << oneNormOfError << "\n"
                 << "||A X - L^-1 AOrig L^-H X||_oo = " 
                 << infNormOfError << "\n"
                 << "||A X - L^-1 AOrig L^-H X||_F  = "
                 << frobNormOfError << endl;
        }
    }
    else
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // triu(B)^-H AOrig triu(B)^-1
        Trsm( LEFT, UPPER, NORMAL, diag, F(1), B, Y );
        Hemm( LEFT, UPPER, F(1), AOrig, Y, F(0), Z );
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), B, Z );
        Hemm( LEFT, UPPER, F(-1), A, X, F(1), Z );
        R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
        R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
        R infNormOfA = HermitianNorm( uplo, A, INFINITY_NORM );
        R frobNormOfA = HermitianNorm( uplo, A, FROBENIUS_NORM );
        R oneNormOfError = Norm( Z, ONE_NORM );
        R infNormOfError = Norm( Z, INFINITY_NORM );
        R frobNormOfError = Norm( Z, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            cout << "||AOrig||_1 = ||AOrig||_oo     = "
                 << infNormOfAOrig << "\n"
                 << "||AOrig||_F                    = "
                 << frobNormOfAOrig << "\n"
                 << "||A||_1 = ||A||_oo             = "
                 << infNormOfA << "\n"
                 << "||A||_F                        = "
                 << frobNormOfA << "\n"
                 << "||A X - U^-H AOrig U^-1 X||_1  = "
                 << oneNormOfError << "\n"
                 << "||A X - U^-H AOrig U^-1 X||_oo = " 
                 << infNormOfError << "\n"
                 << "||A X - U^-H AOrig U^-1 X||_F  = "
                 << frobNormOfError << endl;
        }
    }
}

template<typename F> 
void TestTwoSidedTrsm
( bool testCorrectness, bool print, UpperOrLower uplo, UnitOrNonUnit diag,
  int m, const Grid& g )
{
    DistMatrix<F> A(g), B(g), AOrig(g);

    Zeros( m, m, A );
    Zeros( m, m, B );
    MakeHermitianUniformSpectrum( A, 1, 10 );
    MakeHermitianUniformSpectrum( B, 1, 10 );
    MakeTriangular( uplo, B );
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
    {
        A.Print("A");
        B.Print("B");
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting reduction to Hermitian standard EVP...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TwoSidedTrsm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::val )
        gFlops *= 4.;
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = "
             << gFlops << endl;
    }
    if( print )
        A.Print("A after reduction");
    if( testCorrectness )
        TestCorrectness( print, uplo, diag, A, B, AOrig );
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
        int r = Input("--r","height of process grid",0);
        const char uploChar = Input
            ("--uplo","lower or upper triangular storage: L/U",'L');
        const char diagChar = Input("--unit","(non-)unit diagonal: N/U",'N');
        const int m = Input("--m","height of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test TwoSidedTrsm" << uploChar << diagChar << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestTwoSidedTrsm<double>( testCorrectness, print, uplo, diag, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestTwoSidedTrsm<Complex<double> >
        ( testCorrectness, print, uplo, diag, m, g );
    }
    catch( ArgException& e ) { }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << endl; 
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }
    Finalize();
    return 0;
}
