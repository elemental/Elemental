/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HermitianNorm.hpp"
#include "elemental/lapack-like/HermitianTridiag.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
using namespace std;
using namespace elem;

template<typename R>
void TestCorrectness
( bool print,
  UpperOrLower uplo, 
  const DistMatrix<R>& A, 
        DistMatrix<R>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,STAR> d(g);
    DistMatrix<R,MD,STAR> e(g);
    A.GetDiagonal( d );
    A.GetDiagonal( e, subdiagonal );

    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    // The unaligned [MD,STAR] <- [MD,STAR] redistribution is not yet written,
    // so go around it via [MD,STAR] <- [STAR,STAR] <- [MD,STAR]
    DistMatrix<R,STAR,STAR> e_STAR_STAR(g);
    DistMatrix<R,MD,STAR> eOpposite(g);
    e_STAR_STAR = e;
    eOpposite.AlignWithDiagonal( A, -subdiagonal );
    eOpposite = e_STAR_STAR;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<R> B(g);
    B.AlignWith( A );
    Zeros( m, m, B );
    B.SetDiagonal( d );
    B.SetDiagonal( e, subdiagonal );
    B.SetDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    if( uplo == LOWER )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, B );
        ApplyPackedReflectors
        ( RIGHT, LOWER, VERTICAL, BACKWARD, subdiagonal, A, B );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, subdiagonal, A, B );
        ApplyPackedReflectors
        ( RIGHT, UPPER, VERTICAL, FORWARD, subdiagonal, A, B );
    }

    // Compare the appropriate triangle of AOrig and B
    MakeTriangular( uplo, AOrig );
    MakeTriangular( uplo, B );
    Axpy( R(-1), AOrig, B );

    const R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    const R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    const R infNormOfError = HermitianNorm( uplo, B, INFINITY_NORM );
    const R frobNormOfError = HermitianNorm( uplo, B, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||A - Q^H T Q||_oo         = " << infNormOfError << "\n"
             << "    ||A - Q^H T Q||_F          = " << frobNormOfError << endl;
    }
}

template<typename R> 
void TestCorrectness
( bool print,
  UpperOrLower uplo, 
  const DistMatrix<Complex<R> >& A, 
  const DistMatrix<Complex<R>,STAR,STAR>& t,
        DistMatrix<Complex<R> >& AOrig )
{
    typedef Complex<R> C;
    const Grid& g = A.Grid();
    const int m = AOrig.Height();

    int subdiagonal = ( uplo==LOWER ? -1 : +1 );

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    DistMatrix<R,MD,STAR> d(g);
    DistMatrix<R,MD,STAR> e(g);
    A.GetRealPartOfDiagonal( d );
    A.GetRealPartOfDiagonal( e, subdiagonal );
     
    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    DistMatrix<R,STAR,STAR> e_STAR_STAR(g);
    DistMatrix<R,MD,STAR> eOpposite(g);
    e_STAR_STAR = e;
    eOpposite.AlignWithDiagonal( A.DistData(), -subdiagonal );
    eOpposite = e_STAR_STAR;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<C> B(g);
    B.AlignWith( A );
    Zeros( m, m, B );
    B.SetRealPartOfDiagonal( d );
    B.SetRealPartOfDiagonal( e, subdiagonal );
    B.SetRealPartOfDiagonal( eOpposite, -subdiagonal );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    if( uplo == LOWER )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, 
          UNCONJUGATED, subdiagonal, A, t, B );
        ApplyPackedReflectors
        ( RIGHT, LOWER, VERTICAL, BACKWARD, 
          CONJUGATED, subdiagonal, A, t, B );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, UPPER, VERTICAL, FORWARD, 
          UNCONJUGATED, subdiagonal, A, t, B );
        ApplyPackedReflectors
        ( RIGHT, UPPER, VERTICAL, FORWARD, 
          CONJUGATED, subdiagonal, A, t, B );
    }

    // Compare the appropriate triangle of AOrig and B
    MakeTriangular( uplo, AOrig );
    MakeTriangular( uplo, B );
    Axpy( C(-1), AOrig, B );

    const R infNormOfAOrig = HermitianNorm( uplo, AOrig, INFINITY_NORM );
    const R frobNormOfAOrig = HermitianNorm( uplo, AOrig, FROBENIUS_NORM );
    const R infNormOfError = HermitianNorm( uplo, B, INFINITY_NORM );
    const R frobNormOfError = HermitianNorm( uplo, B, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||AOrig||_1 = ||AOrig||_oo = " << infNormOfAOrig << "\n"
             << "    ||AOrig||_F                = " << frobNormOfAOrig << "\n"
             << "    ||AOrig - Q^H A Q||_oo     = " << infNormOfError << "\n"
             << "    ||AOrig - Q^H A Q||_F      = " << frobNormOfError << endl;
    }
}

template<typename R>
void TestRealHermitianTridiag
( bool testCorrectness, bool print,
  UpperOrLower uplo, int m, const Grid& g )
{
    DistMatrix<R> A(g), AOrig(g);

    HermitianUniformSpectrum( m, A, -10, 10 );
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
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    HermitianTridiag( uplo, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double gFlops = 4./3.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        A.Print("A after HermitianTridiag");
    if( testCorrectness )
        TestCorrectness( print, uplo, A, AOrig );
}

template<typename R>
void TestComplexHermitianTridiag
( bool testCorrectness, bool print,
  UpperOrLower uplo, int m, const Grid& g )
{
    typedef Complex<R> C;
    DistMatrix<C> A(g), AOrig(g);
    DistMatrix<C,STAR,STAR> t(g);

    HermitianUniformSpectrum( m, A, -10, 10 );
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
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    HermitianTridiag( uplo, A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double gFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        A.Print("A after HermitianTridiag");
        t.Print("t after HermitianTridiag");
    }
    if( testCorrectness )
        TestCorrectness( print, uplo, A, t, AOrig );
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
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocal );
        SetLocalSymvBlocksize<Complex<double> >( nbLocal );
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test HermitianTridiag" << uploChar << endl;

        if( commRank == 0 )
        {
            cout << "----------------------------------\n"
                 << "Double-precision normal algorithm:\n"
                 << "----------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestRealHermitianTridiag<double>
        ( testCorrectness, print, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, row-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestRealHermitianTridiag<double>( testCorrectness, print, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------------------\n"
                 << "Double-precision square algorithm, col-major grid:\n"
                 << "--------------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestRealHermitianTridiag<double>( testCorrectness, print, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "------------------------------------------\n"
                 << "Double-precision complex normal algorithm:\n"
                 << "------------------------------------------" << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_NORMAL );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, print, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "row-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( ROW_MAJOR );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, print, uplo, m, g );

        if( commRank == 0 )
        {
            cout << "-------------------------------------------\n"
                 << "Double-precision complex square algorithm, \n"
                 << "col-major grid:\n"
                 << "-------------------------------------------" 
                 << endl;
        }
        SetHermitianTridiagApproach( HERMITIAN_TRIDIAG_SQUARE );
        SetHermitianTridiagGridOrder( COLUMN_MAJOR );
        TestComplexHermitianTridiag<double>
        ( testCorrectness, print, uplo, m, g );
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
