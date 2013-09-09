/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Circulant.hpp"
#include "elemental/matrices/Fourier.hpp"
#include "elemental/matrices/Zeros.hpp"
#include "elemental/io.hpp"
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrices?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create a circulant matrix
        const Grid& g = DefaultGrid();
        std::vector<C> a( n );
        for( Int j=0; j<n; ++j )
            a[j] = j;
        auto A = Circulant( g, a );
        if( display )
            Display( A, "Circulant" );
        if( print )
            Print( A, "Circulant matrix:" );

        // Create a Fourier matrix, which can be used to diagonalize circulant
        // matrices
        auto F = Fourier<Real>( g, n );
        if( display )
            Display( F, "DFT matrix" );
        if( print )
            Print( F, "DFT matrix:" );
        
        // Form B := A F
        auto B = Zeros<C>( g, n, n );
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
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
