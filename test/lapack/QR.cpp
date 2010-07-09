/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Generates random matrix then solves for its QR factorization."
         << endl << endl;
    cout << "  QR <r> <c> <m> <n> <nb> <test correctness?> "
         << "<print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  m: height of matrix            " << endl;
    cout << "  n: width of matrix             " << endl;
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
( bool printMatrices,
  const DistMatrix<R,MC,MR>& A,
        DistMatrix<R,MC,MR>& ARef )
{
    const Grid& g = A.GetGrid();
    const int m = A.Height();
    const int n = A.Width();

    if( g.VCRank() == 0 )
        cout << "  Testing orthogonality of Q:" << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<R,MC,MR> Z(m,n,g);
    Z.SetToIdentity();
    lapack::UT( Left, Lower, ConjugateTranspose, 0, A, Z );
    lapack::UT( Left, Lower, Normal, 0, A, Z );

    DistMatrix<R,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, n, n );

    // Form Identity
    DistMatrix<R,MC,MR> X(n,n,g);
    X.SetToIdentity();

    // Form X := I - Q^H Q
    blas::Axpy( (R)-1, ZUpper, X );

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

    if( g.VCRank() == 0 )
    {
        cout << "  Testing if A = QR...";
        cout.flush();
    }

    // Form Q R
    DistMatrix<R,MC,MR> U( A );
    U.MakeTrapezoidal( Left, Upper );
    lapack::UT( Left, Lower, ConjugateTranspose, 0, A, U );

    // Form Q R - A
    blas::Axpy( (R)-1, ARef, U );
    
    // Compute the maximum deviance
    R myMaxDevFromA = 0.;
    for( int j=0; j<U.LocalWidth(); ++j )
        for( int i=0; i<U.LocalHeight(); ++i )
            myMaxDevFromA = 
                max(myMaxDevFromA,abs(U.LocalEntry(i,j)));
    R maxDevFromA;
    Reduce
    ( &myMaxDevFromA, &maxDevFromA, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation from A is " << maxDevFromA << endl;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void TestCorrectness
( bool printMatrices,
  const DistMatrix<complex<R>,MC,MR  >& A,
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& ARef )
{
    typedef complex<R> C;

    const Grid& g = A.GetGrid();
    const int m = A.Height();
    const int n = A.Width();

    if( g.VCRank() == 0 )
    {
        cout << "  Testing orthogonality of Q...";
        cout.flush();
    }

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<C,MC,MR> Z(m,n,g);
    Z.SetToIdentity();
    lapack::UT( Left, Lower, ConjugateTranspose, 0, A, t, Z );
    lapack::UT( Left, Lower, Normal, 0, A, t, Z );
    
    DistMatrix<C,MC,MR> ZUpper(g);
    ZUpper.View( Z, 0, 0, n, n );

    // Form Identity
    DistMatrix<C,MC,MR> X(m,n,g);
    X.SetToIdentity();

    // Form X := I - Q^H Q
    blas::Axpy( (C)-1, Z, X );

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

    if( g.VCRank() == 0 )
    {
        cout << "  Testing if A = QR...";
        cout.flush();
    }

    // Form Q R
    DistMatrix<C,MC,MR> U( A );
    U.MakeTrapezoidal( Left, Upper );
    lapack::UT( Left, Lower, ConjugateTranspose, 0, A, t, U );

    // Form Q R - A
    blas::Axpy( (C)-1, ARef, U );
    
    // Compute the maximum deviance
    R myMaxDevFromA = 0.;
    for( int j=0; j<U.LocalWidth(); ++j )
        for( int i=0; i<U.LocalHeight(); ++i )
            myMaxDevFromA = 
                max(myMaxDevFromA,abs(U.LocalEntry(i,j)));
    R maxDevFromA;
    Reduce
    ( &myMaxDevFromA, &maxDevFromA, 1, MPI_MAX, 0, g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "max deviation from A is " << maxDevFromA << endl;
}
#endif // WITHOUT_COMPLEX

template<typename T>
void TestQR
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g );

template<>
void TestQR<double>
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g);
    DistMatrix<R,MC,MR> ARef(g);

    A.ResizeTo( m, n );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::QR( A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::QRGFlops<R>( m, n, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, ARef );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestQR< complex<double> >
( bool testCorrectness, bool printMatrices,
  int m, int n, const Grid& g )
{
    typedef complex<double> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR  > A(g);
    DistMatrix<C,MD,Star> t(g);
    DistMatrix<C,MC,MR  > ARef(g);

    A.ResizeTo( m, n );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( g.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( g.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( g.VCRank() == 0 )
    {
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::QR( A, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::QRGFlops<C>( m, n, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( printMatrices, A, t, ARef );
}
#endif

int main( int argc, char* argv[] )
{
    int rank;
    elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 8 )
    {
        if( rank == 0 )
            Usage();
        elemental::Finalize();
        return 0;
    }
    try
    {
        const int   r = atoi( argv[1] );
        const int   c = atoi( argv[2] );
        const int   m = atoi( argv[3] );
        const int   n = atoi( argv[4] );
        const int   nb = atoi( argv[5] );
        const bool  testCorrectness = atoi( argv[6] );
        const bool  printMatrices = atoi( argv[7] );
#ifndef RELEASE
        if( rank == 0 )
        {
            cout << "==========================================" << endl;
            cout << " In debug mode! Performance will be poor! " << endl;
            cout << "==========================================" << endl;
        }
#endif
        if( n > m )
            throw logic_error( "QR only supported when height >= width." );
        Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test QR" << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestQR<double>
        ( testCorrectness, printMatrices, m, n, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestQR<dcomplex>
        ( testCorrectness, printMatrices, m, n, g );
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

