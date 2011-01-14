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
    double* A;
#ifndef WITHOUT_COMPLEX
    DComplex* B;
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
    n = 5;
    colAlignment = 0;
    rowAlignment = 0;
    localHeight = ElementalLocalLength( m, MCRank, colAlignment, r );
    localWidth = ElementalLocalLength( n, MRRank, rowAlignment, c );
    ldim = localHeight;

    A = (double*)malloc(ldim*localWidth*sizeof(double));
#ifndef WITHOUT_COMPLEX
    B = (DComplex*)malloc(ldim*localWidth*sizeof(DComplex));
#endif /* WITHOUT_COMPLEX */

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
#endif /* WITHOUT_COMPLEX */
            }
        }
    }

    /* Attempt to solve a real Hermitian EVP */
    {
        MC_MR_Double AHandle = ElementalRegister_MC_MR_Double
            ( m, n, colAlignment, rowAlignment, A, ldim, g );
        Star_VR_Double wHandle = ElementalCreateEmpty_Star_VR_Double( g );
        MC_MR_Double ZHandle = ElementalCreateEmpty_MC_MR_Double( g );
        ElementalPrint_MC_MR_Double( "A", AHandle );
#ifndef WITHOUT_PMRRR
        ElementalHermitianEigDouble( 'L', AHandle, wHandle, ZHandle, 1 );
        ElementalPrint_Star_VR_Double( "eigenvalues", wHandle );
        ElementalPrint_MC_MR_Double( "eigenvectors", ZHandle );
#endif /* WITHOUT_PMRRR */
    }

#ifndef WITHOUT_COMPLEX
    /* Attempt to solve a complex Hermitian EVP */
    {
        MC_MR_DComplex BHandle = ElementalRegister_MC_MR_DComplex
            ( m, n, colAlignment, rowAlignment, B, ldim, g );
        Star_VR_Double wHandle = ElementalCreateEmpty_Star_VR_Double( g );
        MC_MR_DComplex ZHandle = ElementalCreateEmpty_MC_MR_DComplex( g );
        ElementalPrint_MC_MR_DComplex( "B", BHandle );
#ifndef WITHOUT_PMRRR
        ElementalHermitianEigDComplex( 'L', BHandle, wHandle, ZHandle, 1 );
        ElementalPrint_Star_VR_Double( "eigenvalues", wHandle );
        ElementalPrint_MC_MR_DComplex( "eigenvectors", ZHandle );
#endif /* WITHOUT_PMRRR */
    }
#endif /* WITHOUT_COMPLEX */

    free( A );
#ifndef WITHOUT_COMPLEX
    free( B );
#endif /* WITHOUT_COMPLEX */
    ElementalFinalize();
    return 0;
}

