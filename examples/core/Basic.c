/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-C.h"

int 
main( int argc, char* argv[] )
{
    ElInitialize( &argc, &argv );
    const ElInt m = ElInput_I("--m","matrix height",10);
    const ElInt n = ElInput_I("--n","matrix width",10);
    const ElInt mSub = ElInput_I("--mSub","submatrix height",5);
    const ElInt nSub = ElInput_I("--nSub","submatrix width",5);
    const bool print = ElInput_b("--print","print matrix?",false);
    const bool display = ElInput_b("--display","display matrix?",true);
    ElProcessInput();
    ElPrintInputReport();

    if( mSub > m || nSub > n )
    {
        int worldRank;
        MPI_Comm_rank( MPI_COMM_WORLD, &worldRank );
        if( worldRank == 0 )
            printf("Invalid submatrix dimensions\n");
        ElFinalize();
        return 0;
    }

    const ElGrid* grid = ElGridCreate( MPI_COMM_WORLD, EL_COLUMN_MAJOR );
    ElDistMatrix_d* A = ElDistMatrixCreateSpecific_d( EL_MR, EL_MC, grid );
    ElDistMatrixResize_d( A, m, n );
    
    ElInt i, j;
    for( j=0; j<n; ++j )
        for( i=0; i<m; ++i )
            ElDistMatrixSet_d( A, i, j, i+j );

    if( print )
        ElPrintDistMatrix_d( A, "A" );
    if( display )
        ElDisplayDistMatrix_d( A, "A" );

    /* Extract an mSub x nSub submatrix */
    ElInt* rowInds = malloc(mSub*sizeof(ElInt));
    ElInt* colInds = malloc(nSub*sizeof(ElInt));
    for( i=0; i<mSub; ++i )
        rowInds[i] = rand() % m;
    for( j=0; j<nSub; ++j )
        colInds[j] = rand() % n;
    ElDistMatrix_d* ASub = 
        ElDistMatrixGetSubmatrix_d( A, mSub, rowInds, nSub, colInds );
    if( print )
        ElPrintDistMatrix_d( ASub, "ASub" );
    if( display )
        ElDisplayDistMatrix_d( ASub, "ASub" );

    ElFinalize();
    return 0;
}
