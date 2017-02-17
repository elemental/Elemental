/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.h"

complex_double MapFunc( complex_double alpha )
{ return alpha*alpha; }

int
main( int argc, char* argv[] )
{
    ElError error = ElInitialize( &argc, &argv );
    if( error != EL_SUCCESS )
        MPI_Abort( MPI_COMM_WORLD, 1 );

    ElInt m, n, mSub, nSub;
    bool print, display;
    error = ElInput_I("--m","matrix height",100,&m);
    EL_ABORT_ON_ERROR( error );
    error = ElInput_I("--n","matrix width",100,&n);
    EL_ABORT_ON_ERROR( error );
    error = ElInput_b("--print","print matrix?",true,&print);
    EL_ABORT_ON_ERROR( error );
    error = ElInput_b("--display","display matrix?",false,&display);
    error = ElProcessInput();
    EL_ABORT_ON_ERROR( error );
    error = ElPrintInputReport();
    EL_ABORT_ON_ERROR( error );

    ElGrid grid;
    error = ElGridCreate( MPI_COMM_WORLD, EL_COLUMN_MAJOR, &grid );
    EL_ABORT_ON_ERROR( error );;

    ElDistMatrix_z A;
    error = ElDistMatrixCreateSpecific_z( EL_MR, EL_MC, grid, &A );
    EL_ABORT_ON_ERROR( error );
    error = ElDistMatrixResize_z( A, m, n );
    EL_ABORT_ON_ERROR( error );

    ElInt i, j;
    for( j=0; j<n; ++j )
    {
        for( i=0; i<m; ++i )
        {
            error = ElDistMatrixSet_z( A, i, j, i+j );
            EL_ABORT_ON_ERROR( error );
        }
    }
    if( print )
    {
        error = ElPrintDist_z( A, "A" );
        EL_ABORT_ON_ERROR( error );
    }
    if( display )
    {
        error = ElDisplayDist_z( A, "A" );
        EL_ABORT_ON_ERROR( error );
    }

    ElEntrywiseMapDist_z( A, &MapFunc );
    if( print )
    {
        error = ElPrintDist_z( A, "ASquared" );
        EL_ABORT_ON_ERROR( error );
    }
    if( display )
    {
        error = ElDisplayDist_z( A, "ASquared" );
        EL_ABORT_ON_ERROR( error );
    }

    error = ElFinalize();
    if( error != EL_SUCCESS )
        MPI_Abort( MPI_COMM_WORLD, 1 );
    return 0;
}
