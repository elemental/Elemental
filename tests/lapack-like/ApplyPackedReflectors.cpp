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

template<typename R> 
void TestCorrectness
( LeftOrRight side, UpperOrLower uplo, ForwardOrBackward order,
  int offset, bool printMatrices, const DistMatrix<R>& H )
{
    const Grid& g = H.Grid();
    const int m = H.Height();

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of transform..." << endl;

    // Form Z := Q^H Q or Q^H Q as an approximation to identity
    DistMatrix<R> Y(m,m,g);
    MakeIdentity( Y );
    ApplyPackedReflectors( side, uplo, VERTICAL, order, offset, H, Y );
    if( printMatrices )
    {
        DistMatrix<R> W(m,m,g);
        MakeIdentity( W );
        if( order == FORWARD )
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, BACKWARD, offset, H, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, FORWARD, offset, H, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<R> Z(m,m,g);
    Zero( Z );
    Syrk( uplo, NORMAL, 1.0, Y, 0.0, Z );

    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<R> X(m,m,g);
    MakeIdentity( X );
    Axpy( R(-1), Z, X );
    if( printMatrices )
    {
        if( order == FORWARD )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    const R oneNormOfError = Norm( X, ONE_NORM );
    const R infNormOfError = Norm( X, INFINITY_NORM );
    const R frobNormOfError = Norm( X, FROBENIUS_NORM );
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

template<typename R> 
void TestCorrectness
( LeftOrRight side, UpperOrLower uplo, ForwardOrBackward order,
  Conjugation conjugation, int offset, bool printMatrices,
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t )
{
    typedef Complex<R> C;

    const Grid& g = H.Grid();
    const int m = H.Height();

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of transform..." << endl;

    // Form Z := Q^H Q or Q Q^H as an approximation to identity
    DistMatrix<C> Y(m,m,g);
    MakeIdentity( Y );
    ApplyPackedReflectors
    ( side, uplo, VERTICAL, order, conjugation, offset, H, t, Y );
    if( printMatrices )
    {
        DistMatrix<C> W(m,m,g);
        MakeIdentity( W );
        if( order == FORWARD )
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, BACKWARD, conjugation, offset, H, t, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            ApplyPackedReflectors
            ( side, uplo, VERTICAL, FORWARD, conjugation, offset, H, t, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<C> Z(m,m,g);
    MakeZeros( Z );
    Herk( uplo, NORMAL, C(1), Y, C(0), Z );
    
    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<C> X(m,m,g);
    MakeIdentity( X );
    Axpy( C(-1), Z, X );
    if( printMatrices )
    {
        if( order == FORWARD )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    // Compute the maximum deviance
    const R oneNormOfError = Norm( X, ONE_NORM );
    const R infNormOfError = Norm( X, INFINITY_NORM );
    const R frobNormOfError = Norm( X, FROBENIUS_NORM );
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

template<typename R>
void TestRealUT
( LeftOrRight side, UpperOrLower uplo, 
  ForwardOrBackward order, Conjugation conjugation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    DistMatrix<R> H(g), A(g);
    Uniform( m, m, H );
    Uniform( m, m, A );
    if( printMatrices )
    {
        H.Print("H");
        A.Print("A");
    }

    if( g.Rank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    ApplyPackedReflectors( side, uplo, VERTICAL, order, offset, H, A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double gFlops = 2.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( side, uplo, order, offset, printMatrices, H );
}

template<typename R>
void TestComplexUT
( LeftOrRight side, UpperOrLower uplo, 
  ForwardOrBackward order, Conjugation conjugation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    typedef Complex<R> C;

    DistMatrix<C> H(g), A(g);
    Uniform( m, m, H );
    Uniform( m, m, A );

    const int diagLength = DiagonalLength(H.Height(),H.Width(),offset);
    DistMatrix<C,MD,STAR> t(g);
    t.AlignWithDiagonal( H, offset );
    t.ResizeTo( diagLength, 1 );

    DistMatrix<C> HCol(g);
    if( uplo == LOWER )
    {
        for( int i=0; i<t.Height(); ++i )
        {
            // View below the diagonal containing the implicit 1
            View( HCol, H, i-offset+1, i, m-(i-offset+1), 1 );
            C norm = Nrm2( HCol );
            C alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }
    else
    {
        for( int i=0; i<t.Height(); ++i ) 
        {
            // View above the diagonal containing the implicit 1
            View( HCol, H, 0, i+offset, i, 1 );
            C norm = Nrm2( HCol );
            C alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }

    if( printMatrices )
    {
        H.Print("H");
        A.Print("A");
        t.Print("t");
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
    const double gFlops = 8.*Pow(double(m),3.)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
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
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--gridHeight","height of process grid",0);
        const char sideChar = Input("--side","side to apply from: L/R",'L');
        const char uploChar = Input("--uplo","store in triangle: L/U",'L');
        const bool forward = Input("--forward","forward application?",true);
        const bool conjugate = Input("--conjugate","conjugate?",false);
        const int m = Input("--height","height of matrix",100);
        const int offset = Input("--offset","diagonal offset for storage",0);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness  = Input
            ("--correctness","test correctness?",true);
        const bool printMatrices = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        const LeftOrRight side = CharToLeftOrRight( sideChar );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const ForwardOrBackward order = ( forward ? FORWARD : BACKWARD );
        const Conjugation conjugation = 
            ( conjugate ? CONJUGATED : UNCONJUGATED );
        SetBlocksize( nb );
        if( uplo == LOWER && offset > 0 )
            throw logic_error
            ("Offset cannot be positive if transforms are in lower triangle");
        else if( uplo == UPPER && offset < 0 )
            throw logic_error
            ("Offset cannot be negative if transforms are in upper triangle");
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test UT transform" << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestRealUT<double>
        ( side, uplo, order, conjugation, m, offset, 
          testCorrectness, printMatrices, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestComplexUT<double>
        ( side, uplo, order, conjugation, m, offset, 
          testCorrectness, printMatrices, g );
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
