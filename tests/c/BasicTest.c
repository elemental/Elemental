#include "elemental/exports/c.h"
#include <stdio.h>
#include <stdlib.h>

int 
main( int argc, char* argv[] )
{
    Grid g;
    int r, c;
    int VCRank, MCRank, MRRank;
    int m, n;
    int colAlignment, rowAlignment;
    int localHeight, localWidth, ldim;
    double* ABuffer;
    MC_MR_Double A;
#ifndef WITHOUT_COMPLEX
    DComplex* BBuffer;
    MC_MR_DComplex B;
#endif /* WITHOUT_COMPLEX */

    ElementalInit( &argc, &argv );

    g = CreateDefaultGrid( ELEMENTAL_COMM_WORLD );

    r = GridHeight( g );
    c = GridWidth( g );
    VCRank = GridVCRank( g );
    MCRank = GridMCRank( g );
    MRRank = GridMRRank( g );

    if( VCRank == 0 )
        fprintf( stdout, "Grid is %d x %d.\n", r, c );

    m = 6;
    n = 6;
    colAlignment = 0;
    rowAlignment = 0;
    localHeight = LocalLength( m, MCRank, colAlignment, r );
    localWidth = LocalLength( n, MRRank, rowAlignment, c );
    ldim = localHeight;

    ABuffer = (double*)malloc(ldim*localWidth*sizeof(double));
    A = Register_MC_MR_Double
        ( m, n, colAlignment, rowAlignment, ABuffer, ldim, g );
#ifndef WITHOUT_COMPLEX
    BBuffer = (DComplex*)malloc(ldim*localWidth*sizeof(DComplex));
    B = Register_MC_MR_DComplex
        ( m, n, colAlignment, rowAlignment, BBuffer, ldim, g );
#endif /* WITHOUT_COMPLEX */

    /* Generate matrices that may be interpreted as Hermitian by 
     * (conjugate-)transposing the data from the lower/upper triangle.
     * This only imposes a constraint on complex matrices, where the 
     * diagonal must be real-valued.
     */
    {
        int iLocal, jLocal;
        for( jLocal=0; jLocal<localWidth; ++jLocal )
        {
            int j = MRRank + jLocal*c;
            for( iLocal=0; iLocal<localHeight; ++iLocal )
            {
                int i = MCRank + iLocal*r;
                /* We can actually easily generate a Hermitian matrix in 
                 * parallel by setting entry (i,j) to (i+j,i-j)
                 */
                ABuffer[iLocal+jLocal*ldim] = i+j;
#ifndef WITHOUT_COMPLEX
                BBuffer[iLocal+jLocal*ldim].real = i+j; 
                BBuffer[iLocal+jLocal*ldim].imag = i-j;
#endif /* WITHOUT_COMPLEX */
            }
        }
    }

    /* Print A and B */
    Print_MC_MR_Double( "A", A );
    Print_MC_MR_DComplex( "B", B );

#ifndef WITHOUT_PMRRR
    /* Attempt to solve a real Hermitian EVP */
    {
        Star_VR_Double w = CreateEmpty_Star_VR_Double( g );
        MC_MR_Double Z = CreateEmpty_MC_MR_Double( g );
        HermitianEigDouble( 'L', 1, 'A', 0, 0, 0, 0, A, w, Z );
        Print_Star_VR_Double( "eigenvalues of A", w );
        Print_MC_MR_Double( "eigenvectors of A", Z );
    }
#ifndef WITHOUT_COMPLEX
    /* Attempt to solve a complex Hermitian EVP */
    {
        Star_VR_Double w = CreateEmpty_Star_VR_Double( g );
        MC_MR_DComplex Z = CreateEmpty_MC_MR_DComplex( g );
        HermitianEigDComplex( 'L', 1, 'A', 0, 0, 0, 0, B, w, Z );
        Print_Star_VR_Double( "eigenvalues of B", w );
        Print_MC_MR_DComplex( "eigenvectors of B", Z );
    }
#endif /* WITHOUT_COMPLEX */
#endif /* WITHOUT_PMRRR */

    free( ABuffer );
#ifndef WITHOUT_COMPLEX
    free( BBuffer );
#endif /* WITHOUT_COMPLEX */
    ElementalFinalize();
    return 0;
}

