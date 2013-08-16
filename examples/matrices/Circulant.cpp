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
        DistMatrix<Complex<double>> A;
        std::vector<Complex<double>> a( n );
        for( Int j=0; j<n; ++j )
            a[j] = j;
        Circulant( A, a );
        if( display )
            Display( A, "Circulant" );
        if( print )
            Print( A, "Circulant matrix:" );

        // Create a Fourier matrix, which can be used to diagonalize circulant
        // matrices
        DistMatrix<Complex<double>> F;
        Fourier( F, n );
        if( display )
            Display( F, "DFT matrix" );
        if( print )
            Print( F, "DFT matrix:" );
        
        // Form B := A F
        DistMatrix<Complex<double>> B;
        Zeros( B, n, n );
        Gemm( NORMAL, NORMAL, Complex<double>(1), A, F, Complex<double>(0), B );

        // Form A := F^H B = F^H \hat A F
        Gemm
        ( ADJOINT, NORMAL, Complex<double>(1), F, B, Complex<double>(0), A );
        if( display )
            Display( A, "F^H A F" );
        if( print )
            Print( A, "A := F^H A F" );

        // Form the thresholded result
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        for( Int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( Int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                const double absValue = Abs(A.GetLocal(iLocal,jLocal));
                if( absValue < 1e-13 )
                    A.SetLocal(iLocal,jLocal,0);
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
