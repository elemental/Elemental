/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Tests UT transform application." << endl << endl;
    cout << "  UT <r> <c> <side> <shape> <orientation> <m> <offset> <nb> "
         << "<correctness?> <print?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  side: {L/R}                    " << endl;
    cout << "  shape: {L/U}                   " << endl;
    cout << "  orientation: {N/C}             " << endl;
    cout << "  m: height of matrix            " << endl;
    cout << "  offset: diagonal transforms are stored above/below" << endl;
    cout << "  nb: algorithmic blocksize      " << endl;
    cout << "  test correctness?: false iff 0 " << endl;
    cout << "  print matrices?: false iff 0   " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( double truth, double computed )
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-10 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-10 ); }
#endif

template<typename R>
void TestCorrectness
( Side side, 
  Shape shape,
  Orientation orientation,
  int offset,
  bool printMatrices,
  const DistMatrix<R,MC,MR>& H )
{
    const Grid& g = H.GetGrid();
    const int m = H.Height();

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of transform:" << endl;

    // Form Z := Q^H Q or Q^H Q as an approximation to identity
    DistMatrix<R,MC,MR> Y(m,m,g);
    Y.SetToIdentity();
    lapack::UT( side, shape, orientation, offset, H, Y );
    if( printMatrices )
    {
        DistMatrix<R,MC,MR> W(m,m,g);
        W.SetToIdentity();
        if( orientation == Normal )
        {
            lapack::UT( side, shape, ConjugateTranspose, offset, H, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            lapack::UT( side, shape, Normal, offset, H, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<R,MC,MR> Z(m,m,g);
    blas::Syrk( Lower, Normal, 1.0, Y, 0.0, Z );

    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<R,MC,MR> X(m,m,g);
    X.SetToIdentity();
    blas::Axpy( (R)-1, Z, X );
    if( printMatrices )
    {
        if( orientation == Normal )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    // Compute the maximum deviance
    R myMaxDevFromIdentity = 0.;
    for( int j=0; j<X.LocalWidth(); ++j )
        for( int i=0; i<X.LocalHeight(); ++i )
            myMaxDevFromIdentity = 
                max(myMaxDevFromIdentity,abs(X.LocalEntry(i,j)));
    R maxDevFromIdentity;
    Reduce
    ( &myMaxDevFromIdentity, &maxDevFromIdentity, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation from I is " << maxDevFromIdentity << endl;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void TestCorrectness
( Side side,
  Shape shape,
  Orientation orientation,
  int offset,
  bool printMatrices,
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t )
{
    typedef complex<R> C;

    const Grid& g = H.GetGrid();
    const int m = H.Height();

    if( g.VCRank() == 0 )
    {
        cout << "  Testing orthogonality of transform...";
        cout.flush();
    }

    // Form Z := Q^H Q or Q Q^H as an approximation to identity
    DistMatrix<C,MC,MR> Y(m,m,g);
    Y.SetToIdentity();
    lapack::UT( side, shape, orientation, offset, H, t, Y );
    if( printMatrices )
    {
        DistMatrix<C,MC,MR> W(m,m,g);
        W.SetToIdentity();
        if( orientation == Normal )
        {
            lapack::UT( side, shape, ConjugateTranspose, offset, H, t, W );
            Y.Print("Q");
            W.Print("Q^H");
        }
        else
        {
            lapack::UT( side, shape, Normal, offset, H, t, W );
            Y.Print("Q^H");
            W.Print("Q");
        }
    }
    DistMatrix<C,MC,MR> Z(m,m,g);
    blas::Herk( Lower, Normal, (C)1, Y, (C)0, Z );
    
    // Form X := I - Q^H Q or Q Q^H
    DistMatrix<C,MC,MR> X(m,m,g);
    X.SetToIdentity();
    blas::Axpy( (C)-1, Z, X );
    if( printMatrices )
    {
        if( orientation == Normal )
            X.Print("I - Q Q^H");
        else
            X.Print("I - Q^H Q");
    }

    // Compute the maximum deviance
    R myMaxDevFromIdentity = 0.;
    for( int j=0; j<X.LocalWidth(); ++j )
        for( int i=0; i<X.LocalHeight(); ++i )
            myMaxDevFromIdentity = 
                max(myMaxDevFromIdentity,abs(X.LocalEntry(i,j)));
    R maxDevFromIdentity;
    Reduce
    ( &myMaxDevFromIdentity, &maxDevFromIdentity, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation from I is " << maxDevFromIdentity << endl;
}
#endif // WITHOUT_COMPLEX

template<typename T>
void TestUT
( Side side, Shape shape, Orientation orientation, 
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g );

template<>
void TestUT<double>
( Side side, Shape shape, Orientation orientation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> H(g);
    DistMatrix<R,MC,MR> A(g);

    H.ResizeTo( m, m );
    A.ResizeTo( m, m );

    H.SetToRandom();
    A.SetToRandom();
    if( printMatrices )
    {
        H.Print("H");
        A.Print("A");
    }

    if( g.VCRank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::UT( side, shape, orientation, offset, H, A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::UTGFlops<R>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( side, shape, orientation, offset, printMatrices, H );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestUT< complex<double> >
( Side side, Shape shape, Orientation orientation,
  int m, int offset, bool testCorrectness, bool printMatrices,
  const Grid& g )
{
    typedef complex<double> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR  > H(g);
    DistMatrix<C,MD,Star> t(g);
    DistMatrix<C,MC,MR  > A(g);

    H.ResizeTo( m, m );
    A.ResizeTo( m, m );

    t.AlignWithDiag( H, offset );
    t.ResizeTo( H.DiagonalLength( offset ), 1 );

    H.SetToRandom();
    A.SetToRandom();
    DistMatrix<C,MC,MR> HCol(g);
    if( shape == Lower )
    {
        for( int i=0; i<t.Height(); ++i )
        {
            // View below the diagonal containing the implicit 1
            HCol.View( H, i-offset+1, i, m-(i-offset+1), 1 );
            C norm = blas::Nrm2( HCol );
            C alpha = 2./(norm*norm+1.);
            t.Set( i, 0, alpha );
        }
    }
    else
    {
        for( int i=0; i<t.Height(); ++i ) 
        {
            // View above the diagonal containing the implicit 1
            HCol.View( H, 0, i+offset, i, 1 );
            C norm = blas::Nrm2( HCol );
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

    if( g.VCRank() == 0 )
    {
        cout << "  Starting UT transform...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::UT( side, shape, orientation, offset, H, t, A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::UTGFlops<C>( m, runTime );
    if( g.VCRank() == 0 )
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
        ( side, shape, orientation, offset, printMatrices, H, t );
    }
}
#endif

int 
main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 11 )
    {
        if( rank == 0 )
            Usage();
        elemental::Finalize();
        return 0;
    }
    try
    {
        const int         r = atoi( argv[1] );
        const int         c = atoi( argv[2] );
        const Side        side = CharToSide( *argv[3] );
        const Shape       shape = CharToShape( *argv[4] );
        const Orientation orientation = CharToOrientation( *argv[5] );
        const int         m = atoi( argv[6] );
        const int         offset = atoi( argv[7] );
        const int         nb = atoi( argv[8] );
        const bool        testCorrectness = atoi( argv[9] );
        const bool        printMatrices = atoi( argv[10] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test UT transform" << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestUT<double>
        ( side, shape, orientation, m, offset, 
          testCorrectness, printMatrices, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestUT<dcomplex>
        ( side, shape, orientation, m, offset, 
          testCorrectness, printMatrices, g );
        if( rank == 0 )
            cout << endl;
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:" << endl 
             << e.what() << endl;
    }   
    elemental::Finalize();
    return 0;
}

