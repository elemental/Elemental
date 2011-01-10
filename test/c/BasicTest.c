#include "elemental/exports/c.h"
#include <stdio.h>

int main( int argc, char* argv[] )
{
    ElementalInit( &argc, &argv );

    int gridHandle = ElementalDefaultGrid( MPI_COMM_WORLD );

    int r = ElementalGridHeight( gridHandle );
    int c = ElementalGridWidth( gridHandle );
    int VCRank = ElementalGridVCRank( gridHandle );

    fprintf
    ( stdout, "Process %d: Grid is %d x %d.\n", VCRank, r, c );

    ElementalFinalize();
    return 0;
}

