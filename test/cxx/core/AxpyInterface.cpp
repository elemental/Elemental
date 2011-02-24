#include "elemental.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::imports;

void Usage()
{
    cout << "AxpyInterface <n>\n"
         << "\n"
         << "  <n>: size of matrix to test.\n"
         << std::endl;
}

// TODO: Make this more interesting
void
FormDistributedMatrix( DistMatrix<double,MC,MR>& A )
{
    const int m = A.Height();
    const int n = A.Width();
    Matrix<double> X( m, n );
    X.SetToRandom();

    AxpyInterface<double> interface( A );
    if( A.Grid().VCRank() == 0 )
        interface.Axpy( 1.0, X, 0, 0 );
}

int
main( int argc, char* argv[] )
{
    Init( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );

    if( argc < 2 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try 
    {
        const int n = atoi(argv[1]);

        Grid g( comm );

        DistMatrix<double,MC,MR> A( n, n, g );
        A.SetToRandom();
        FormDistributedMatrix( A );
        A.Print("A");

        DistMatrix<int,VC,Star> p_VC_Star( n, 1, g );
        DistMatrix<double,MC,MR> LUOfA( A );
        advanced::LU( LUOfA, p_VC_Star );
        LUOfA.Print("LU(A)");

        DistMatrix<int,Star,Star> p_Star_Star( g );
        p_Star_Star = p_VC_Star;
        vector<int> image, preimage;
        advanced::internal::ComposePivots( p_Star_Star, image, preimage, 0 );
        advanced::internal::ApplyRowPivots( A, image, preimage, 0 );
        A.Print("A after pivoting");

        DistMatrix<double,MC,MR> b( n, 1, g );
        b.SetToRandom();
        b.Print("b");

        basic::Trsm
        ( Left, Lower, Normal, Unit, 1.0, LUOfA, b );
        basic::Trsm
        ( Left, Upper, Normal, NonUnit, 1.0, LUOfA, b );
        b.Print("x");

        DistMatrix<double,MC,MR> z( n, 1, g );
        basic::Gemm( Normal, Normal, 1.0, A, b, 0.0, z );
        z.Print("PAx ?= b");
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n"
             << e.what() << endl;
    }

    Finalize();
    return 0;
}

