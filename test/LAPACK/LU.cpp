#include <cmath>
#include <ctime>
#include <sstream>
#include "Elemental.h"
#include "ElementalLAPACKInternal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

void Usage()
{
    cout << "Generates random matrix then solves for its LU factors."
         << endl << endl;
    cout << "  LU <r> <c> <m> <nb> <test correctness?> "
         << "<print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  m: height of matrix            " << endl;
    cout << "  nb: algorithmic blocksize      " << endl;
    cout << "  test correctness?: false iff 0 " << endl;
    cout << "  print matrices?: false iff 0   " << endl;
    cout << endl;
}

template<typename T>
bool OKRelativeError( T truth, T computed );

template<>
bool OKRelativeError( float truth, float computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(float)1) <= 1e-3 );
}

template<>
bool OKRelativeError( double truth, double computed )
{
    return ( fabs(truth-computed) / max(fabs(truth),(double)1) <= 1e-10 );
}

#ifndef WITHOUT_COMPLEX
template<>
bool OKRelativeError( scomplex truth, scomplex computed )
{
    return ( norm(truth-computed) / max(norm(truth),(float)1) <= 1e-3 );
}

template<>
bool OKRelativeError( dcomplex truth, dcomplex computed )
{
    return ( norm(truth-computed) / max(norm(truth),(double)1) <= 1e-10 );
}
#endif

template<typename T>
void TestCorrectness
(       DistMatrix<T,Star,Star>& A_ref,
  const DistMatrix<T,MC,MR>& A,
  const DistMatrix<int,VC,Star>& p,
  const bool printMatrices             )
{
    const Grid& grid = A.GetGrid();
    const int m = A_ref.Height();
    DistMatrix<T,Star,Star> A_copy(grid);
    DistMatrix<int,Star,Star> p_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    A_copy = A;
    p_copy = p;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    DistMatrix<int,Star,Star> p_ref(m,1,grid);
    LAPACK::LU( A_ref.LocalMatrix(), p_ref.LocalMatrix() );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
    {
        A_ref.Print("True A:");
        p_ref.Print("True p:");
    }

    if( grid.VCRank() == 0 )
    {
        cout << "  Testing correctness...";
        cout.flush();
    }
    // Test A
    for( int j=0; j<m; ++j )
    {
        for( int i=0; i<m; ++i )
        {
            T truth = A_ref.LocalEntry(i,j);
            T computed = A_copy.LocalEntry(i,j);

            if( ! OKRelativeError( truth, computed ) )
            {
                ostringstream msg;
                msg << "FAILED at index (" << i << "," << j << "): truth=" 
                     << truth << ", computed=" << computed;
                throw msg.str();
            }
        }
    }
    for( int i=0; i<m; ++i )
    {
        int truth = p_ref.LocalEntry(i,0);
        int computed = p_copy.LocalEntry(i,0);
        if( truth != computed+1 /* 0 vs. 1 indexing */ )
        {
            ostringstream msg;
            msg << "Pivots off at index " << i << ": truth=" << truth 
                 << ", computed=" << computed;
            throw msg.str();
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestLU
( const int m, const bool testCorrectness, const bool printMatrices, 
  const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,Star,Star> A_ref(grid);
    DistMatrix<int,VC,Star> p(grid);

    A.ResizeTo( m, m );
    p.ResizeTo( m, 1 );

    A.SetToRandom();
    if( testCorrectness )
    {
        if( grid.VCRank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        A_ref = A;
        if( grid.VCRank() == 0 )
            cout << "DONE" << endl;
    }
    if( printMatrices )
    {
        A.Print("A");
    }

    if( grid.VCRank() == 0 )
    {
        cout << "  Starting LU factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    LAPACK::LU( A, p );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = LAPACK::Internal::LUGFlops<T>( m, runTime );
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        A.Print("A after factorization");
        p.Print("p after factorization");
    }
    if( testCorrectness )
    {
        TestCorrectness( A_ref, A, p, printMatrices );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 7 )
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
        const int   m = atoi( argv[3] );
        const int   nb = atoi( argv[4] );
        const bool  testCorrectness = atoi( argv[5] );
        const bool  printMatrices = atoi( argv[6] );
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
            cout << "Will test LU" << endl;

        if( m<=200 )
        {
            if( rank == 0 )
            {
                cout << "--------------------" << endl;
                cout << "Testing with floats:" << endl;
                cout << "--------------------" << endl;
            }
            TestLU<float>
            ( m, testCorrectness, printMatrices, grid );
            if( rank == 0 )
                cout << endl;
        }
        else
        {
            if( rank == 0 )
            {
                cout << "--------------------------------" << endl;
                cout << "Floats unsuitable for this test." << endl;
                cout << "--------------------------------" << endl;
                cout << endl;
            }
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        TestLU<double>
        ( m, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;

#ifndef WITHOUT_COMPLEX
        if( m <= 200 )
        {
            if( rank == 0 )
            {
                cout << "--------------------------------------" << endl;
                cout << "Testing with single-precision complex:" << endl;
                cout << "--------------------------------------" << endl;
            }
            TestLU<scomplex>
            ( m, testCorrectness, printMatrices, grid );
            if( rank == 0 )
                cout << endl;
        }
        else
        {
            if( rank == 0 )
            {
                cout << "----------------------------------------" << endl;
                cout << "Complex floats unsuitable for this test." << endl;
                cout << "----------------------------------------" << endl;
                cout << endl;
            }
        }

        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        TestLU<dcomplex>
        ( m, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;
#endif
    }
    catch( string errorMsg )
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

