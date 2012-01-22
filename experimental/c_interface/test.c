#include "./elemental.h"
#include <stdio.h>
#include <stdlib.h>

int
main( int argc, char* argv[] )
{    
    /* Handles for Elemental C++ objects */
    GridHandle grid;
    RealDistMatHandle A, B, X;
    RealDistColVecHandle w;

    /* Process grid information */
    int gridHeight, gridWidth, gridRow, gridCol, gridRank;

    /* Dimensions of our local matrix (for A and B) */
    int localHeight, localWidth; 

    /* Local buffers for the distributed A and B matrices */
    double *ABuffer, *BBuffer;

    /* Indices */
    int i, j, iLocal, jLocal;

    /* Useful constants */
    int n=10;                       /* problem size */
    int nb=96;                      /* algorithmic blocksize */
    MPI_Comm comm = MPI_COMM_WORLD; /* global communicator */

    /* Initialize Elemental and MPI */
    Initialize( &argc, &argv );

    /* Create a process grid and extract the relevant information */
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
    for( jLocal=0; jLocal<localWidth; ++jLocal )
    {
        j = gridCol + jLocal*gridWidth;
        for( iLocal=0; iLocal<localHeight; ++iLocal )
        {
            i = gridRow + iLocal*gridHeight;
            ABuffer[iLocal+jLocal*localHeight] = (double)i+j;
        }
    }

    /* Set B to twice the identity since it is a trivial SPD matrix */
    for( jLocal=0; jLocal<localWidth; ++jLocal )
    {
        j = gridCol + jLocal*gridWidth;
        for( iLocal=0; iLocal<localHeight; ++iLocal )
        {
            i = gridRow + iLocal*gridHeight;
            if( i == j )
                BBuffer[iLocal+jLocal*localHeight] = 2.0;
            else
                BBuffer[iLocal+jLocal*localHeight] = 0.0;
        }
    }

    /* Register the distributed matrices, A and B, with Elemental */
    A = RegisterRealDistMat( n, n, 0, 0, ABuffer, localHeight, grid );
    B = RegisterRealDistMat( n, n, 0, 0, BBuffer, localHeight, grid );

    /* Print the input matrices */
    if( gridRank == 0 )
        printf("A:\n");
    PrintRealDistMat( A );
    if( gridRank == 0 )
        printf("B:\n");
    PrintRealDistMat( B );

    /* Set the algorithmic blocksize to 'nb' */
    SetBlocksize( nb );

    /* Run the eigensolver */
    if( gridRank == 0 )
        printf("Solving for (w,X) in AX=BXW\n");
    SymmetricAxBx( A, B, &w, &X );

    /* Print the eigenvalues and eigenvectors */
    if( gridRank == 0 )
        printf("X:\n");
    PrintRealDistMat( X );
    if( gridRank == 0 )
        printf("w:\n");
    PrintRealDistColVec( w );

    /* Shut down Elemental and MPI */
    Finalize();

    return 0;
}
