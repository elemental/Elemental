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
    cout << "Generates SPD matrix then solves for its Cholesky factor."
         << endl << endl;
    cout << "  Chol <r> <c> <shape> <var2/3> <m> <nb> <test correctness?> "
         << "<print matrices?>" << endl << endl;
    cout << "  r: number of process rows      " << endl;
    cout << "  c: number of process cols      " << endl;
    cout << "  shape: {L,U}                   " << endl;
    cout << "  var2/3: 2 iff 0                " << endl;
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
( const Shape shape,
        DistMatrix<T,Star,Star>& A_ref,
  const DistMatrix<T,MC,MR>& A,
  const bool printMatrices                       )
{
    const Grid& grid = A.GetGrid();
    const int m = A_ref.Height();
    DistMatrix<T,Star,Star> A_copy(grid);

    if( grid.VCRank() == 0 )
    {
        cout << "  Gathering computed result...";
        cout.flush();
    }
    A_copy = A;
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( grid.VCRank() == 0 )
    {
        cout << "  Computing 'truth'...";
        cout.flush();
    }
    LAPACK::Chol( shape, A_ref.LocalMatrix() );
    if( grid.VCRank() == 0 )
        cout << "DONE" << endl;

    if( printMatrices )
        A_ref.Print("Truth");

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
                T truth = A_ref.LocalEntry(i,j);
                T computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j << "): truth=" 
                         << truth << ", computed=" << computed;
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
                T truth = A_ref.LocalEntry(i,j);
                T computed = A_copy.LocalEntry(i,j);

                if( ! OKRelativeError( truth, computed ) )
                {
                    ostringstream msg;
                    msg << "FAILED at index (" << i << "," << j << "): truth=" 
                         << truth << ", computed=" << computed;
                    const string s = msg.str();
                    throw s.c_str();
                }
            }
        }
    }
    Barrier( grid.VCComm() );
    if( grid.VCRank() == 0 )
        cout << "PASSED" << endl;
}

template<typename T>
void TestChol
( const Shape shape, const bool var3, const int m, 
  const bool testCorrectness, const bool printMatrices, const Grid& grid )
{
    double startTime, endTime, runTime, gFlops;
    DistMatrix<T,MC,MR> A(grid);
    DistMatrix<T,Star,Star> A_ref(grid);

    A.ResizeTo( m, m );

    A.SetToRandomDiagDominant();
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
        cout << "  Starting Cholesky factorization...";
        cout.flush();
    }
    Barrier( MPI_COMM_WORLD );
    startTime = Time();
    if( var3 )
        LAPACK::Internal::CholVar3( shape, A );
    else
        LAPACK::Internal::CholVar2( shape, A );
    Barrier( MPI_COMM_WORLD );
    endTime = Time();
    runTime = endTime - startTime;
    gFlops = LAPACK::Internal::CholGFlops<T>( m, runTime );
    if( grid.VCRank() == 0 )
        cout << "DONE. GFlops = " << gFlops << endl;
    if( printMatrices )
    {
        A.Print("A after factorization");
    }
    if( testCorrectness )
    {
        TestCorrectness( shape, A_ref, A, printMatrices );
    }
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 9 )
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
        const bool  var3 = ( atoi(argv[4]) != 0 );
        const int   m = atoi( argv[5] );
        const int   nb = atoi( argv[6] );
        const bool  testCorrectness = atoi( argv[7] );
        const bool  printMatrices = atoi( argv[8] );
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
            cout << "Will test Chol" << ShapeToChar(shape) << ", Var"
                 << ( var3 ? "3" : "2" ) << endl;

        if( m<=200 )
        {
            if( rank == 0 )
            {
                cout << "--------------------" << endl;
                cout << "Testing with floats:" << endl;
                cout << "--------------------" << endl;
            }
            TestChol<float>
            ( shape, var3, m, testCorrectness, printMatrices, grid );
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
        TestChol<double>
        ( shape, var3, m, testCorrectness, printMatrices, grid );
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
            TestChol<scomplex>
            ( shape, var3, m, testCorrectness, printMatrices, grid );
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
        TestChol<dcomplex>
        ( shape, var3, m, testCorrectness, printMatrices, grid );
        if( rank == 0 )
            cout << endl;
#endif
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

