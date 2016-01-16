/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrices?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create a circulant matrix
        vector<C> a( n );
        for( Int j=0; j<n; ++j )
            a[j] = j;
        DistMatrix<C> A;
        Circulant( A, a );
        if( display )
            Display( A, "Circulant" );
        if( print )
            Print( A, "Circulant matrix:" );

        // Create a Fourier matrix, which can be used to diagonalize circulant
        // matrices
        DistMatrix<C> F;
        Fourier( F, n );
        if( display )
            Display( F, "DFT matrix" );
        if( print )
            Print( F, "DFT matrix:" );
        
        // Form B := A F
        DistMatrix<C> B;
        Zeros( B, n, n );
        Gemm( NORMAL, NORMAL, C(1), A, F, C(0), B );

        // Form A := F^H B = F^H \hat A F
        Gemm( ADJOINT, NORMAL, C(1), F, B, C(0), A );
        if( display )
            Display( A, "F^H A F" );
        if( print )
            Print( A, "A := F^H A F" );

        // Form the thresholded result
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const double absValue = Abs(A.GetLocal(iLoc,jLoc));
                if( absValue < 1e-13 )
                    A.SetLocal(iLoc,jLoc,0);
            }
        }
        if( display )
            Display( A, "Thresholded (1e-13) A" );
        if( print )
            Print( A, "A with values below 1e-13 removed" );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
