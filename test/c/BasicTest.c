#include "elemental/exports/c.h"
#include <stdio.h>
#include <stdlib.h>

int main( int argc, char* argv[] )
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
#endif // WITHOUT_COMPLEX

    ElementalInit( &argc, &argv );

    g = ElementalDefaultGrid( MPI_COMM_WORLD );

    r = ElementalGridHeight( g );
    c = ElementalGridWidth( g );
    VCRank = ElementalGridVCRank( g );
    MCRank = ElementalGridMCRank( g );
    MRRank = ElementalGridMRRank( g );

    if( VCRank == 0 )
        fprintf( stdout, "Grid is %d x %d.\n", r, c );

    m = 6;
    n = 6;
    colAlignment = 0;
    rowAlignment = 0;
    localHeight = ElementalLocalLength( m, MCRank, colAlignment, r );
    localWidth = ElementalLocalLength( n, MRRank, rowAlignment, c );
    ldim = localHeight;

    ABuffer = (double*)malloc(ldim*localWidth*sizeof(double));
    A = ElementalRegister_MC_MR_Double
        ( m, n, colAlignment, rowAlignment, ABuffer, ldim, g );
#ifndef WITHOUT_COMPLEX
    BBuffer = (DComplex*)malloc(ldim*localWidth*sizeof(DComplex));
    B = ElementalRegister_MC_MR_DComplex
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

    /* Attempt to solve a real Hermitian EVP */
    {
        Star_VR_Double w = ElementalCreateEmpty_Star_VR_Double( g );
        MC_MR_Double Z = ElementalCreateEmpty_MC_MR_Double( g );
        ElementalPrint_MC_MR_Double( "A", A );
#ifndef WITHOUT_PMRRR
        ElementalHermitianEigDouble( 'L', A, w, Z, 1 );
        ElementalPrint_Star_VR_Double( "eigenvalues", w );
        ElementalPrint_MC_MR_Double( "eigenvectors", Z );
#endif /* WITHOUT_PMRRR */
    }

#ifndef WITHOUT_COMPLEX
    /* Attempt to solve a complex Hermitian EVP */
    {
        Star_VR_Double w = ElementalCreateEmpty_Star_VR_Double( g );
        MC_MR_DComplex Z = ElementalCreateEmpty_MC_MR_DComplex( g );
        ElementalPrint_MC_MR_DComplex( "B", B );
#ifndef WITHOUT_PMRRR
        ElementalHermitianEigDComplex( 'L', B, w, Z, 1 );
        ElementalPrint_Star_VR_Double( "eigenvalues of (lower) Hermitian B", w );
        ElementalPrint_MC_MR_DComplex( "eigenvectors of (lower) Hermitian B", Z );
#endif /* WITHOUT_PMRRR */
    }
#endif /* WITHOUT_COMPLEX */

    free( ABuffer );
#ifndef WITHOUT_COMPLEX
    free( BBuffer );
#endif /* WITHOUT_COMPLEX */
    ElementalFinalize();
    return 0;
}

