#include "./c_interface.hpp"
#include <stdlib.h>

int
main( int argc, char* argv[] )
{
    int gridHeight, gridWidth, gridRow, gridCol, gridRank;
    int n=10, localHeight, localWidth, i, j, iLocal, jLocal;
    double *ABuffer, *BBuffer;
    MPI_Comm comm;
    GridHandle grid;
    RealDistMatHandle A, B, X;
    RealDistColVecHandle w;

    Initialize( &argc, &argv );
    comm = MPI_COMM_WORLD;
    grid = CreateGrid( comm );
    gridHeight = GridHeight( grid );
    gridWidth = GridWidth( grid );
    gridRow = GridRow( grid );
    gridCol = GridCol( grid );
    gridRank = GridRank( grid );

    /* Create buffers for passing into data for distributed matrices */
    localHeight = LocalLength( n, gridRow, gridHeight );
    localWidth = LocalLength( n, gridCol, gridWidth );
    ABuffer = (double*)malloc(localHeight*localWidth*sizeof(double));
    BBuffer = (double*)malloc(localHeight*localWidth*sizeof(double));

    /* Set entry (i,j) of the A matrix to i+j, which is symmetric */
    for( jLocal=gridCol; jLocal<n; jLocal+=gridWidth )
    {
        j = gridCol + jLocal*gridWidth;
        for( iLocal=gridRow; iLocal<n; iLocal+=gridHeight )
        {
            i = gridRow + iLocal*gridHeight;
            ABuffer[iLocal+jLocal*localHeight] = (double)i+j;
        }
    }

    /* Set B to twice the identity since it is a trivial SPD matrix */
    for( jLocal=gridCol; jLocal<n; jLocal+=gridWidth )
    {
        j = gridCol + jLocal*gridWidth;
        for( iLocal=gridRow; iLocal<n; iLocal+=gridHeight )
        {
            i = gridRow + iLocal*gridHeight;
            if( i == j )
                ABuffer[iLocal+jLocal*localHeight] = 2.0;
            else
                ABuffer[iLocal+jLocal*localHeight] = 0.0;
        }
    }

    /* Register the distributed matrices, A and B, with Elemental */
    A = RegisterRealDistMat( n, n, 0, 0, ABuffer, localHeight, grid );
    B = RegisterRealDistMat( n, n, 0, 0, BBuffer, localHeight, grid );

    if( gridRank == 0 )
        printf("A:\n");
    PrintRealDistMat( A );
    if( gridRank == 0 )
        printf("B:\n");
    PrintRealDistMat( B );

    /* Create empty objects for the eigenvalues, w, and eigenvectors, X */
    w = CreateEmptyRealDistColVec( grid );
    X = CreateEmptyRealDistMat( grid );

    if( gridRank == 0 )
        printf("Solving for (w,X) in AX=BXW\n");
    SymmetricAxBx( A, B, w, X );

    if( gridRank == 0 )
        printf("X:");
    PrintRealDistMat( X );
    if( gridRank == 0 )
        printf("w:");
    PrintRealDistColVec( w );

    Finalize();
    return 0;
}
