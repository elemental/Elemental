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
#ifndef WITHOUT_COMPLEX
    ElementalDComplex* B;
#endif // WITHOUT_COMPLEX

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
#ifndef WITHOUT_COMPLEX
    B = (ElementalDComplex*)malloc(ldim*localWidth*sizeof(ElementalDComplex));
#endif // WITHOUT_COMPLEX

    /* Mark the entries owned by each process */
    {
        int iLocal, jLocal;
        for( jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( iLocal=0; iLocal<localHeight; ++iLocal )
            {
                A[iLocal+jLocal*ldim] = VCRank;
#ifndef WITHOUT_COMPLEX
                B[iLocal+jLocal*ldim].real = MCRank;
                B[iLocal+jLocal*ldim].imag = MRRank;
#endif // WITHOUT_COMPLEX
            }
        }
    }

    /* Print the real distributed matrix */
    {
        int AHandle = ElementalDistMatrixDouble
            ( m, n, colAlignment, rowAlignment, A, ldim, gridHandle );
        ElementalDistMatrixDoublePrint( "A", AHandle );
    }

#ifndef WITHOUT_COMPLEX
    /* Print the complex distributed matrix */
    {
        int BHandle = 
            ElementalDistMatrixDComplex
            ( m, n, colAlignment, rowAlignment, B, ldim, gridHandle );
        ElementalDistMatrixDComplexPrint( "B", BHandle );
    }
#endif // WITHOUT_COMPLEX

    free( A );
#ifndef WITHOUT_COMPLEX
    free( B );
#endif // WITHOUT_COMPLEX
    ElementalFinalize();
    return 0;
}

