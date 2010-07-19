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
{ return ( Abs(truth-computed) / max(Abs(truth),(double)1) <= 1e-8 ); }

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{ return ( Abs(truth-computed) / max(Abs(truth),(double)1) <= 1e-8 ); }
#endif

template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<R,MC,MR>& A, 
        DistMatrix<R,Star,Star>& ARef )
{
    const Grid& g = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<R,Star,Star> ACopy(g);

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    ACopy = A;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    double startTime = Time();
    lapack::Tridiag( shape, ARef.LocalMatrix() );
    double stopTime = Time();
    double gFlops = lapack::internal::TridiagGFlops<R>(m,stopTime-startTime);
    if( g.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;

    if( printMatrices )
        ARef.Print("True A:");

    if( g.VCRank() == 0 )
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
                R computed = ACopy.LocalEntry(i,j);

                if( !OKRelativeError(truth,computed) )
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
                R computed = ACopy.LocalEntry(i,j);

                if( !OKRelativeError(truth,computed) )
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
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<complex<R>,MC,  MR  >& A, 
  const DistMatrix<complex<R>,MD,  Star>& t,
        DistMatrix<complex<R>,Star,Star>& ARef )
{
    typedef complex<R> C;
    const Grid& g = A.GetGrid();
    const int m = ARef.Height();
    Matrix<C> tRef;
    DistMatrix<C,Star,Star> ACopy(g);
    

    if( g.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    ACopy = A;
    if( g.VCRank() == 0 )
        cout << "DONE" << endl;

    if( g.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    double startTime = Time();
    lapack::Tridiag( shape, ARef.LocalMatrix(), tRef );
    double stopTime = Time();
    double gFlops = lapack::internal::TridiagGFlops<C>(m,stopTime-startTime);
    if( g.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;

    if( printMatrices )
        ARef.Print("True A:");

    if( g.VCRank() == 0 )
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
                C computed = ACopy.LocalEntry(i,j);

                if( !OKRelativeError(truth,computed) )
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
                C computed = ACopy.LocalEntry(i,j);

                if( !OKRelativeError(truth,computed) )
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
    Barrier( g.VCComm() );
    if( g.VCRank() == 0 )
        cout << "PASSED" << endl;
}
#endif // WITHOUT_COMPLEX

template<typename T>
void TestTridiag
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g );

template<>
void TestTridiag<double>
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<R,MC,MR> A(g);
    DistMatrix<R,Star,Star> ARef(g);

    A.ResizeTo( m, m );

    A.SetToRandomHPD();
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
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops<R>( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
        A.Print("A after Tridiag");
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, ARef );
}

#ifndef WITHOUT_COMPLEX
template<>
void TestTridiag< complex<double> >
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& g )
{
    typedef double R;
    typedef complex<R> C;

    double startTime, endTime, runTime, gFlops;
    DistMatrix<C,MC,MR> A(g);
    DistMatrix<C,MD,Star> t(g);
    DistMatrix<C,Star,Star> ARef(g);

    A.ResizeTo( m, m );

    // Make A diagonally dominant
    A.SetToRandomHPD();
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
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    lapack::Tridiag( shape, A, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = lapack::internal::TridiagGFlops< complex<R> >( m, runTime );
    if( g.VCRank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( printMatrices )
    {
        A.Print("A after Tridiag");
        t.Print("t after Tridiag");
    }
    if( testCorrectness )
        TestCorrectness( printMatrices, shape, A, t, ARef );
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
        Grid g( MPI_COMM_WORLD, r, c );
        SetBlocksize( nb );

        if( rank == 0 )
            cout << "Will test Tridiag" << ShapeToChar(shape) << endl;

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestTridiag<double>( testCorrectness, printMatrices, shape, m, g );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "----------------------------" << endl;
            cout << "Testing with double-complex:" << endl;
            cout << "----------------------------" << endl;
        }
        TestTridiag<dcomplex>( testCorrectness, printMatrices, shape, m, g );
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

