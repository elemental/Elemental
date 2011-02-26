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
        A.SetToZero();

        AxpyInterface<double> interface;
        interface.Attach( A );
        Matrix<double> X( n, n );
        const int rank = A.Grid().VCRank();
        for( int j=0; j<n; ++j )
            for( int i=0; i<n; ++i )
                X.Set(i,j,rank);
        interface.Axpy( 1.0, X, 0, 0 );
        interface.Detach();

        // Ensure that our local matrix is the sum of all the ranks
        const int p = mpi::CommSize( mpi::COMM_WORLD );
        const double sumOfRanks = ((p-1)*(p-1)+(p-1))/2;
        // Check that our local matrix is equal to sumOfRanks everywhere
        double myMaxError = 0.;
        for( int jLocal=0; jLocal<A.LocalWidth(); ++jLocal )
            for( int iLocal=0; iLocal<A.LocalHeight(); ++iLocal )
                myMaxError = 
                    std::max( myMaxError, 
                              Abs(sumOfRanks-A.GetLocalEntry(iLocal,jLocal)) );
        double maxError; 
        mpi::AllReduce
        ( &myMaxError, &maxError, 1, mpi::SUM, g.VCComm() );

        if( rank == 0 )
            std::cout << "max error = " << maxError << std::endl;
        if( maxError > 0.000001 )
            A.Print("A");
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

