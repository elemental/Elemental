#include "elemental/exports/c.h"
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
{
    int gridHandle;
    int r, c;
    int VCRank, MCRank, MRRank;
    int m, n;
    int colAlignment, rowAlignment;
    int localHeight, localWidth, ldim;
    double* A;
    ElementalDComplex* B;

    ElementalInit( &argc, &argv );

    gridHandle = ElementalDefaultGrid( MPI_COMM_WORLD );

    r = ElementalGridHeight( gridHandle );
    c = ElementalGridWidth( gridHandle );
    VCRank = ElementalGridVCRank( gridHandle );
    MCRank = ElementalGridMCRank( gridHandle );
    MRRank = ElementalGridMRRank( gridHandle );

    if( VCRank == 0 )
        fprintf( stdout, "Grid is %d x %d.\n", r, c );

    m = 6;
    n = 5;
    colAlignment = 0;
    rowAlignment = 0;
    localHeight = ElementalLocalLength( m, MCRank, colAlignment, r );
    localWidth = ElementalLocalLength( n, MRRank, rowAlignment, c );
    ldim = localHeight;

    A = (double*)malloc(ldim*localWidth*sizeof(double));
    B = (ElementalDComplex*)malloc(ldim*localWidth*sizeof(ElementalDComplex));

    /* Mark the entries owned by each process */
    {
        int iLocal, jLocal;
        for( jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( iLocal=0; iLocal<localHeight; ++iLocal )
            {
                A[iLocal+jLocal*ldim] = VCRank;
                B[iLocal+jLocal*ldim].real = MCRank;
                B[iLocal+jLocal*ldim].imag = MRRank;
            }
        }
    }

    /* Print the distributed matrices */
    {
        int AHandle, BHandle;
        AHandle = 
            ElementalDistMatrix_MC_MR_Double
            ( m, n, colAlignment, rowAlignment, A, ldim, gridHandle );
        BHandle = 
            ElementalDistMatrix_MC_MR_DComplex
            ( m, n, colAlignment, rowAlignment, B, ldim, gridHandle );

        ElementalDistMatrixDoublePrint( "A", AHandle );
        ElementalDistMatrixDComplexPrint( "B", BHandle );
    }

    free( A );
    free( B );
    ElementalFinalize();
    return 0;
}

