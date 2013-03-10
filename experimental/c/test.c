/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "./elemental.h"
#include <stdio.h>
#include <stdlib.h>

int
main( int argc, char* argv[] )
{    
    /* Handles for Elemental C++ objects */
    ElemGrid grid;
    ElemDistMat A, B, X;
    ElemDistMat_VR_STAR w;

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
    ElemInitialize( &argc, &argv );

    /* Create a process grid and extract the relevant information */
    grid = ElemCreateGrid( comm );
    gridHeight = ElemGridHeight( grid );
    gridWidth = ElemGridWidth( grid );
    gridRow = ElemGridRow( grid );
    gridCol = ElemGridCol( grid );
    gridRank = ElemGridRank( grid );

    /* Create buffers for passing into data for distributed matrices */
    localHeight = ElemLength( n, gridRow, gridHeight );
    localWidth = ElemLength( n, gridCol, gridWidth );
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
    A = ElemRegisterDistMat( n, n, 0, 0, ABuffer, localHeight, grid );
    B = ElemRegisterDistMat( n, n, 0, 0, BBuffer, localHeight, grid );

    /* Print the input matrices */
    if( gridRank == 0 )
        printf("A:\n");
    ElemPrintDistMat( A );
    if( gridRank == 0 )
        printf("B:\n");
    ElemPrintDistMat( B );

    /* Set the algorithmic blocksize to 'nb' */
    ElemSetBlocksize( nb );

    /* Run the eigensolver */
    if( gridRank == 0 )
        printf("Solving for (w,X) in AX=BXW\n");
    ElemSymmetricAxBx( A, B, &w, &X );

    /* Print the eigenvalues and eigenvectors */
    if( gridRank == 0 )
        printf("X:\n");
    ElemPrintDistMat( X );
    if( gridRank == 0 )
        printf("w:\n");
    ElemPrintDistMat_VR_STAR( w );

    /* Shut down Elemental and MPI */
    ElemFinalize();

    return 0;
}
