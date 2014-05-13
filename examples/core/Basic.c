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
    const ElInt m = ElInput_I( "--m", "matrix height", 10 );
    const ElInt n = ElInput_I( "--n", "matrix width", 10 );
    ElProcessInput();
    ElPrintInputReport();

    ElMatrix_d* A = ElMatrixCreate_d();
    ElMatrixResize_d( A, m, n );
    
    ElInt i, j;
    for( j=0; j<n; ++j )
        for( i=0; i<m; ++i )
            ElMatrixSet_d( A, i, j, i+j );

    ElPrintMatrix_d( A, "A" );

    ElFinalize();
    return 0;
}
