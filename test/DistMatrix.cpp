#include <cstdlib>
#include <ctime>
#include "Elemental.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

void Usage()
{
    cout << "Test all of the aligned redistributions of DistMatrix class."
         << endl << endl;
    cout << "  DistMatrix <r> <c> <m> <n>  " << endl << endl;    
    cout << "  r: number of process rows " << endl;
    cout << "  c: number of process cols " << endl;
    cout << "  m: height of matrices     " << endl;
    cout << "  n: width of matrices      " << endl;
    cout << endl;
}

template<typename T, Distribution AColDist, Distribution ARowDist,
                     Distribution BColDist, Distribution BRowDist>
void
Check( DistMatrix<T,AColDist,ARowDist>& A, 
       DistMatrix<T,BColDist,BRowDist>& B )
{
    const Grid& grid = A.GetGrid();

    int p = grid.Size();
    int rank = grid.VCRank();
    int height = B.Height();
    int width = B.Width();
    DistMatrix<T,Star,Star> A_Star_Star(grid);
    DistMatrix<T,Star,Star> B_Star_Star(grid);

    if( rank == 0 )
    {
        cout << "Testing [" << DistToString(AColDist) << ","
                            << DistToString(ARowDist) << "]"
             << " <- ["     << DistToString(BColDist) << ","
                            << DistToString(BRowDist) << "]...";
        cout.flush();
    }
    A = B;

    A_Star_Star = A;
    B_Star_Star = B;

    bool sameData = true;
    for( int j=0; j<width; ++j )
    {
        for( int i=0; i<height; ++i )
        {
            if( A_Star_Star.LocalEntry(i,j) != B_Star_Star.LocalEntry(i,j) )
            {
                sameData = false;
                break;
            }
        }
        if( ! sameData )
            break;
    }

    bool everyonePassed = true;
    int send = sameData;
    int* recvBuf = new int[p];
    AllGather( &send, 1, recvBuf, 1, grid.VCComm() );
    for( int i=0; i<p; ++i )
    {
        if( ! recvBuf[i] )
        {
            everyonePassed = false;
            break;
        }
    }
    delete[] recvBuf;

    if( rank == 0 )
    {
        if( everyonePassed )
            cout << "PASSED" << endl;
        else
            cout << "FAILED" << endl;
    }
}

template<typename T>
void
DistMatrixTest( int m, int n, const Grid& grid )
{
    DistMatrix<T,MC,  MR  > A_MC_MR(grid);
    DistMatrix<T,MC,  Star> A_MC_Star(grid);
    DistMatrix<T,Star,MR  > A_Star_MR(grid);
    DistMatrix<T,MR,  MC  > A_MR_MC(grid);
    DistMatrix<T,MR,  Star> A_MR_Star(grid);
    DistMatrix<T,Star,MC  > A_Star_MC(grid);
    DistMatrix<T,VC,  Star> A_VC_Star(grid);
    DistMatrix<T,Star,VC  > A_Star_VC(grid);
    DistMatrix<T,VR,  Star> A_VR_Star(grid);
    DistMatrix<T,Star,VR  > A_Star_VR(grid);
    DistMatrix<T,Star,Star> A_Star_Star(grid);

    // Communicate from A[MC,MR] 
    A_MC_MR.ResizeTo( m, n );
    A_MC_MR.SetToRandom();
    Check( A_MC_Star,   A_MC_MR );
    Check( A_Star_MR,   A_MC_MR );
    Check( A_MR_MC,     A_MC_MR );
    Check( A_MR_Star,   A_MC_MR );
    Check( A_Star_MC,   A_MC_MR );
    Check( A_VC_Star,   A_MC_MR );
    Check( A_Star_VC,   A_MC_MR );
    Check( A_VR_Star,   A_MC_MR );
    Check( A_Star_VR,   A_MC_MR );
    Check( A_Star_Star, A_MC_MR );

    // Communicate from A[MC,*]
    A_MC_Star.ResizeTo( m, n );
    A_MC_Star.SetToRandom();
    Check( A_MC_MR,     A_MC_Star );
    Check( A_Star_MR,   A_MC_Star );
    Check( A_MR_MC,     A_MC_Star );
    Check( A_MR_Star,   A_MC_Star );
    Check( A_Star_MC,   A_MC_Star );
    Check( A_VC_Star,   A_MC_Star );
    Check( A_Star_VC,   A_MC_Star );
    Check( A_VR_Star,   A_MC_Star );
    Check( A_Star_VR,   A_MC_Star );
    Check( A_Star_Star, A_MC_Star );

    // Communicate from A[*,MR]
    A_Star_MR.ResizeTo( m, n );
    A_Star_MR.SetToRandom();
    Check( A_MC_MR,     A_Star_MR );
    Check( A_MC_Star,   A_Star_MR );
    Check( A_MR_MC,     A_Star_MR );
    Check( A_MR_Star,   A_Star_MR );
    Check( A_Star_MC,   A_Star_MR );
    Check( A_VC_Star,   A_Star_MR );
    Check( A_Star_VC,   A_Star_MR );
    Check( A_VR_Star,   A_Star_MR );
    Check( A_Star_VR,   A_Star_MR );
    Check( A_Star_Star, A_Star_MR );
    
    // Communicate from A[MR,MC]
    A_MR_MC.ResizeTo( m, n );
    A_MR_MC.SetToRandom();
    Check( A_MC_MR,     A_MR_MC );
    Check( A_MC_Star,   A_MR_MC );
    Check( A_Star_MR,   A_MR_MC );
    Check( A_MR_Star,   A_MR_MC );
    Check( A_Star_MC,   A_MR_MC );
    Check( A_VC_Star,   A_MR_MC );
    Check( A_Star_VC,   A_MR_MC );
    Check( A_VR_Star,   A_MR_MC );
    Check( A_Star_VR,   A_MR_MC );
    Check( A_Star_Star, A_MR_MC );

    // Communicate from A[MR,*]
    A_MR_Star.ResizeTo( m, n );
    A_MR_Star.SetToRandom();
    Check( A_MC_MR,     A_MR_Star );
    Check( A_MC_Star,   A_MR_Star );
    Check( A_Star_MR,   A_MR_Star );
    Check( A_MR_MC,     A_MR_Star );
    Check( A_Star_MC,   A_MR_Star );
    Check( A_VC_Star,   A_MR_Star );
    Check( A_Star_VC,   A_MR_Star );
    Check( A_VR_Star,   A_MR_Star );
    Check( A_Star_VR,   A_MR_Star );
    Check( A_Star_Star, A_MR_Star );

    // Communicate from A[*,MC]
    A_Star_MC.ResizeTo( m, n );
    A_Star_MC.SetToRandom();
    Check( A_MC_MR,     A_Star_MC );
    Check( A_MC_Star,   A_Star_MC );
    Check( A_Star_MR,   A_Star_MC );
    Check( A_MR_MC,     A_Star_MC );
    Check( A_MR_Star,   A_Star_MC );
    Check( A_VC_Star,   A_Star_MC );
    Check( A_Star_VC,   A_Star_MC );
    Check( A_VR_Star,   A_Star_MC );
    Check( A_Star_VR,   A_Star_MC );
    Check( A_Star_Star, A_Star_MC );
 
    // Communicate from A[VC,*]
    A_VC_Star.ResizeTo( m, n );
    A_VC_Star.SetToRandom();
    Check( A_MC_MR,     A_VC_Star );
    Check( A_MC_Star,   A_VC_Star );
    Check( A_Star_MR,   A_VC_Star );
    Check( A_MR_MC,     A_VC_Star );
    Check( A_MR_Star,   A_VC_Star );
    Check( A_Star_MC,   A_VC_Star );
    Check( A_Star_VC,   A_VC_Star );
    Check( A_VR_Star,   A_VC_Star );
    Check( A_Star_VR,   A_VC_Star );
    Check( A_Star_Star, A_VC_Star );

    // Communicate from A[*,VC]
    A_Star_VC.ResizeTo( m, n );
    A_Star_VC.SetToRandom();
    Check( A_MC_MR,     A_Star_VC );
    Check( A_MC_Star,   A_Star_VC );
    Check( A_Star_MR,   A_Star_VC );
    Check( A_MR_MC,     A_Star_VC );
    Check( A_MR_Star,   A_Star_VC );
    Check( A_Star_MC,   A_Star_VC );
    Check( A_VC_Star,   A_Star_VC );
    Check( A_VR_Star,   A_Star_VC );
    Check( A_Star_VR,   A_Star_VC );
    Check( A_Star_Star, A_Star_VC );

    // Communicate from A[VR,*]
    A_VR_Star.ResizeTo( m, n );
    A_VR_Star.SetToRandom();
    Check( A_MC_MR,     A_VR_Star );
    Check( A_MC_Star,   A_VR_Star );
    Check( A_Star_MR,   A_VR_Star );
    Check( A_MR_MC,     A_VR_Star );
    Check( A_MR_Star,   A_VR_Star );
    Check( A_Star_MC,   A_VR_Star );
    Check( A_VC_Star,   A_VR_Star );
    Check( A_Star_VC,   A_VR_Star );
    Check( A_Star_VR,   A_VR_Star );
    Check( A_Star_Star, A_VR_Star );

    // Communicate from A[*,VR]
    A_Star_VR.ResizeTo( m, n );
    A_Star_VR.SetToRandom();
    Check( A_MC_MR,     A_Star_VR );
    Check( A_MC_Star,   A_Star_VR );
    Check( A_Star_MR,   A_Star_VR );
    Check( A_MR_MC,     A_Star_VR );
    Check( A_MR_Star,   A_Star_VR );
    Check( A_Star_MC,   A_Star_VR );
    Check( A_VC_Star,   A_Star_VR );
    Check( A_Star_VC,   A_Star_VR );
    Check( A_VR_Star,   A_Star_VR );
    Check( A_Star_Star, A_Star_VR );

    // Communicate from A[*,*]
    A_Star_Star.ResizeTo( m, n );
    A_Star_Star.SetToRandom();
    Check( A_MC_MR,   A_Star_Star );
    Check( A_MC_Star, A_Star_Star );
    Check( A_Star_MR, A_Star_Star );
    Check( A_MR_MC,   A_Star_Star );
    Check( A_MR_Star, A_Star_Star );
    Check( A_Star_MC, A_Star_Star );
    Check( A_VC_Star, A_Star_Star );
    Check( A_Star_VC, A_Star_Star );
    Check( A_VR_Star, A_Star_Star );
    Check( A_Star_VR, A_Star_Star );
}

int main( int argc, char* argv[] )
{
    int rank;
    Elemental::Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    if( argc != 5 )
    {
        if( rank == 0 )
            Usage();
        Elemental::Finalize();
        return 0;
    }
    int r = atoi(argv[1]);
    int c = atoi(argv[2]);
    int m = atoi(argv[3]);
    int n = atoi(argv[4]);
#ifndef RELEASE
    if( rank == 0 )
    {
        cout << "==========================================" << endl;
        cout << " In debug mode! Performance will be poor! " << endl;
        cout << "==========================================" << endl;
    }
#endif
    try
    {
        const Grid grid( MPI_COMM_WORLD, r, c );

        if( rank == 0 )
        {
            cout << "--------------------" << endl;
            cout << "Testing with floats:" << endl;
            cout << "--------------------" << endl;
        }
        DistMatrixTest<float>( m, n, grid );
        if( rank == 0 )
        {
            cout << endl;
        }

        if( rank == 0 )
        {
            cout << "---------------------" << endl;
            cout << "Testing with doubles:" << endl;
            cout << "---------------------" << endl;
        }
        DistMatrixTest<double>( m, n, grid );
        if( rank == 0 )
        {
            cout << endl;
        }

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with single-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        DistMatrixTest<scomplex>( m, n, grid );
        if( rank == 0 )
        {
            cout << endl;
        }
        
        if( rank == 0 )
        {
            cout << "--------------------------------------" << endl;
            cout << "Testing with double-precision complex:" << endl;
            cout << "--------------------------------------" << endl;
        }
        DistMatrixTest<dcomplex>( m, n, grid );
        if( rank == 0 )
        {
            cout << endl;
        }
#endif
    }
    catch( exception e )
    {
        cerr << "Caught exception on process " << rank << endl;
    }
    Elemental::Finalize();
    return 0;
}
