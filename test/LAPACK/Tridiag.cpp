#include <cmath>
#include <ctime>
#include <sstream>
#include "Elemental.hpp"
#include "Elemental/LAPACKInternal.hpp"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

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
{ return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-9 ); }

template<typename T>
void TestCorrectness
( bool printMatrices,
  Shape shape, 
  const DistMatrix<T,MC,MR>& A, 
  const DistMatrix<T,MD,Star>& d,
  const DistMatrix<T,MD,Star>& e,
  const DistMatrix<T,MD,Star>& t,
        DistMatrix<T,Star,Star>& ARef )
{
    const Grid& grid = A.GetGrid();
    const int m = ARef.Height();
    DistMatrix<T,Star,Star> A_copy(grid);
    DistMatrix<T,Star,Star> d_copy(grid);
    DistMatrix<T,Star,Star> e_copy(grid);
    DistMatrix<T,Star,Star> t_copy(grid);
    DistMatrix<T,Star,Star> dRef(m,1,grid);
    DistMatrix<T,Star,Star> eRef(m-1,1,grid);
    DistMatrix<T,Star,Star> tRef(m-1,1,grid);

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
    if( shape == Lower )
    {
        LAPACK::Tridiag( Lower, ARef.LocalMatrix(), 
                                dRef.LocalMatrix(),
                                eRef.LocalMatrix(),
                                tRef.LocalMatrix() );
    }
    else
    {
        Matrix<T> ATransRef;

        BLAS::Trans( ARef.LockedLocalMatrix(), ATransRef );
        LAPACK::Tridiag( Lower, ATransRef,
                                dRef.LocalMatrix(),
                                eRef.LocalMatrix(),
                                tRef.LocalMatrix() );
        BLAS::Trans( ATransRef, ARef.LocalMatrix() );
    }
    double stopTime = Time();
    double gFlops = LAPACK::Internal::TridiagGFlops<T>(m,stopTime-startTime);
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
                T truth = ARef.LocalEntry(i,j);
                T computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed=" 
                         << computed;
                    const string s = msg.str();
                    throw s.c_str();
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
                T truth = ARef.LocalEntry(i,j);
                T computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j 
                         << ") of A: truth=" << truth << ", computed="
                         << computed;
                    const string s = msg.str();
                    throw s.c_str();
                }
            }
        }
    }
    for( int j=0; j<m; ++j )
    {
        T truth = dRef.LocalEntry(j,0);
        T computed = d_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of d: truth=" << truth
                 << ", computed=" << computed;
            const string s = msg.str();
            throw s.c_str();
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        T truth = eRef.LocalEntry(j,0);
        T computed = e_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of e: truth=" << truth
                 << ", computed=" << computed;
            const string s = msg.str();
            throw s.c_str();
        }
    }
    for( int j=0; j<m-1; ++j )
    {
        T truth = tRef.LocalEntry(j,0);
        T computed = t_copy.LocalEntry(j,0);

        if( ! OKRelativeError( truth, computed ) )
        {
            ostringstream msg;
            msg << "FAILED at index " << j << " of t: truth=" << truth
                 << ", computed=" << computed;
            const string s = msg.str();
            throw s.c_str();
        }
    }

    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestTridiag
( bool testCorrectness, bool printMatrices,
  Shape shape, int m, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,MD,Star> d(grid);
    DistMatrix<T,MD,Star> e(grid);
    DistMatrix<T,MD,Star> t(grid);
    DistMatrix<T,Star,Star> ARef(grid);

    A.ResizeTo( m, m );

    d.AlignWithDiag( A );
    if( shape == Lower )
        e.AlignWithDiag( A, -1 );
    else
        e.AlignWithDiag( A, +1 );
    t.AlignWithDiag( A );

    d.ResizeTo( m,   1 );
    e.ResizeTo( m-1, 1 );
    t.ResizeTo( m-1, 1 );

    // Make A diagonally dominant
    A.SetToRandomDiagDominant();
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
    {
        A.Print("A");
    }

    if( grid.VCRank() == 0 )
    {
        cout << "  Starting tridiagonalization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    LAPACK::Tridiag( shape, A, d, e, t );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = LAPACK::Internal::TridiagGFlops<T>( m, runTime );
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        A.Print("A after Tridiag");
        d.Print("d after Tridiag");
        e.Print("e after Tridiag");
        t.Print("t after Tridiag");
    }
    if( testCorrectness )
    {
        TestCorrectness( printMatrices, shape, A, d, e, t, ARef );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 8 )
    {
        if( rank == 0 )
            Usage();
        Elemental::Finalize();
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
    }
    catch( const char* errorMsg )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:" << endl 
             << errorMsg << endl;
    }   
    Elemental::Finalize();
    return 0;
}

