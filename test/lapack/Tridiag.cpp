/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include <ctime>
#include "elemental.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::wrappers::mpi;

void Usage()
{
    cout << "Tridiagonalizes a symmetric matrix." << endl << endl;
    cout << "  Tridiag <r> <c> <shape> <m> <nb> <test correctness?> "
         << "<print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  shape: {L,U}                   " << endl;
    cout << "  m: height of matrix            " << endl;
    cout << "  nb: algorithmic blocksize      " << endl;
    cout << "  test correctness?: false iff 0 " << endl;
    cout << "  print matrices?: false iff 0   " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( double truth, double computed )
{ return ( Abs(truth-computed) / max(Abs(truth),(double)1) <= 1e-9 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( Abs(truth-computed) / max(Abs(truth),(double)1) <= 1e-9 ); }
#endif

template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<R,MC,MR>& A, 
  const DistMatrix<R,MD,Star>& d,
  const DistMatrix<R,MD,Star>& e,
  const DistMatrix<R,MD,Star>& t,
        DistMatrix<R,Star,Star>& ARef )
{
    const Grid& grid = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<R,Star,Star> A_copy(grid);
    DistMatrix<R,Star,Star> d_copy(grid);
    DistMatrix<R,Star,Star> e_copy(grid);
    DistMatrix<R,Star,Star> t_copy(grid);
    DistMatrix<R,Star,Star> dRef(m,1,grid);
    DistMatrix<R,Star,Star> eRef(m-1,1,grid);
    DistMatrix<R,Star,Star> tRef(m-1,1,grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    A_copy = A;
    d_copy = d;
    e_copy = e;
    t_copy = t;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    double startTime = Time();
    lapack::Tridiag
    ( shape, 
      ARef.LocalMatrix(), 
      dRef.LocalMatrix(), 
      eRef.LocalMatrix(), 
      tRef.LocalMatrix() );
    double stopTime = Time();
    double gFlops = lapack::internal::TridiagGFlops<R>(m,stopTime-startTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;

    if( printMatrices )
    {
        ARef.Print("True A:");
        dRef.Print("True d:");
        eRef.Print("True e:");
        tRef.Print("True t:");
    }

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    if( shape == Lower )
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=j; i<m; ++i )
            {
                R truth = ARef.LocalEntry(i,j);
                R computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed=" 
                         << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    else
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                R truth = ARef.LocalEntry(i,j);
                R computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed="
                         << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    for( int j=0; j<m; ++j )
    {
        R truth = dRef.LocalEntry(j,0);
        R computed = d_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of d: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        R truth = eRef.LocalEntry(j,0);
        R computed = e_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of e: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        R truth = tRef.LocalEntry(j,0);
        R computed = t_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of t: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }

    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<complex<R>,MC,  MR  >& A, 
  const DistMatrix<R,         MD,  Star>& d,
  const DistMatrix<R,         MD,  Star>& e,
  const DistMatrix<complex<R>,MD,  Star>& t,
        DistMatrix<complex<R>,Star,Star>& ARef )
{
    typedef complex<R> C;

    const Grid& grid = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<C,Star,Star> A_copy(grid);
    DistMatrix<R,Star,Star> d_copy(grid);
    DistMatrix<R,Star,Star> e_copy(grid);
    DistMatrix<C,Star,Star> t_copy(grid);
    DistMatrix<R,Star,Star> dRef(m,1,grid);
    DistMatrix<R,Star,Star> eRef(m-1,1,grid);
    DistMatrix<C,Star,Star> tRef(m-1,1,grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    A_copy = A;
    d_copy = d;
    e_copy = e;
    t_copy = t;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    double startTime = Time();
    lapack::Tridiag
    ( shape, 
      ARef.LocalMatrix(), 
      dRef.LocalMatrix(),
      eRef.LocalMatrix(),
      tRef.LocalMatrix() );
    double stopTime = Time();
    double gFlops = lapack::internal::TridiagGFlops<C>(m,stopTime-startTime);
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;

    if( printMatrices )
    {
        ARef.Print("True A:");
        dRef.Print("True d:");
        eRef.Print("True e:");
        tRef.Print("True t:");
    }

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    if( shape == Lower )
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=j; i<m; ++i )
            {
                C truth = ARef.LocalEntry(i,j);
                C computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed=" 
                         << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    else
    {
        for( int j=0; j<m; ++j )
        {
            for( int i=0; i<=j; ++i )
            {
                C truth = ARef.LocalEntry(i,j);
                C computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed="
                         << computed;
                    throw logic_error( msg.str() );
                }
            }
        }
    }
    for( int j=0; j<m; ++j )
    {
        C truth = dRef.LocalEntry(j,0);
        C computed = d_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of d: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        C truth = eRef.LocalEntry(j,0);
        C computed = e_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of e: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        C truth = tRef.LocalEntry(j,0);
        C computed = t_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of t: truth=" << truth
                 << ", computed=" << computed;
            throw logic_error( msg.str() );
        }
    }

    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}
#endif // WITHOUT_COMPLEX

template<typename R>
void TestTridiag
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& grid );

template<>
void TestTridiag<double>
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& grid )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(grid);
    DistMatrix<R,MD,Star> d(grid);
    DistMatrix<R,MD,Star> e(grid);
    DistMatrix<R,MD,Star> t(grid);
    DistMatrix<R,Star,Star> ARef(grid);

    A.ResizeTo( m, m );

    d.AlignWithDiag( A );
    if( shape == Lower )
        e.AlignWithDiag( A, -1 );
    else
        e.AlignWithDiag( A, +1 );

    if( shape == Lower )
    {
        // Lower Tridiag traverses down the diagonal, so we must align t 
        // with the top-left (m-1) x (m-1) submatrix of A.
        t.AlignWithDiag( A );
    }
    else
    {
        // Upper Tridiag traverses up the diagonal, so we must align t 
        // with the bottom-right (m-1) x (m-1) submatrix of A.
        DistMatrix<R,MC,MR> ABR(grid);
        ABR.View( A, 1, 1, m-1, m-1 );
        t.AlignWithDiag( ABR );
    }

    d.ResizeTo( m,   1 );
    e.ResizeTo( m-1, 1 );
    t.ResizeTo( m-1, 1 );

    A.SetToRandomHPD();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( grid.VCRank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A, d, e, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops<R>( m, runTime );
    if( grid.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after Tridiag");
        d.Print("d after Tridiag");
        e.Print("e after Tridiag");
        t.Print("t after Tridiag");
    }
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, d, e, t, ARef );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestTridiag< complex<double> >
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& grid )
{
    typedef double R;
    typedef complex<R> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR> A(grid);
    DistMatrix<R,MD,Star> d(grid);
    DistMatrix<R,MD,Star> e(grid);
    DistMatrix<C,MD,Star> t(grid);
    DistMatrix<C,Star,Star> ARef(grid);

    A.ResizeTo( m, m );

    d.AlignWithDiag( A );
    if( shape == Lower )
        e.AlignWithDiag( A, -1 );
    else
        e.AlignWithDiag( A, +1 );
    if( shape == Lower )
    {
        // Lower Tridiag traverses down the diagonal, so we must align t 
        // with the top-left (m-1) x (m-1) submatrix of A.
        t.AlignWithDiag( A );
    }
    else
    {
        // Upper Tridiag traverses up the diagonal, so we must align t 
        // with the bottom-right (m-1) x (m-1) submatrix of A.
        DistMatrix<C,MC,MR> ABR(grid);
        ABR.View( A, 1, 1, m-1, m-1 );
        t.AlignWithDiag( ABR );
    }

    d.ResizeTo( m,   1 );
    e.ResizeTo( m-1, 1 );
    t.ResizeTo( m-1, 1 );

    // Make A diagonally dominant
    A.SetToRandomHPD();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        ARef = A;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
        A.Print("A");

    if( grid.VCRank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A, d, e, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops< complex<R> >( m, runTime );
    if( grid.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after Tridiag");
        d.Print("d after Tridiag");
        e.Print("e after Tridiag");
        t.Print("t after Tridiag");
    }
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, d, e, t, ARef );
}
#endif // WITHOUT_COMPLEX

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
        const Shape shape = CharToShape( *argv[3] );
        const int   m = atoi( argv[4] );
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
        Grid grid( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test Tridiag" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestTridiag<double>( testCorrectness, printMatrices, shape, m, grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "----------------------------" << endl;
            cout << "Testing with double-complex:" << endl;
            cout << "----------------------------" << endl;
        }
        TestTridiag<dcomplex>( testCorrectness, printMatrices, shape, m, grid );
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

