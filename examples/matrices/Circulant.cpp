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
#include "elemental/graphics.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of matrix",10);
        const bool print = Input("--print","print matrices?",true);
#ifdef HAVE_QT5
        const bool display = Input("--display","display matrices?",true);
#endif
        ProcessInput();
        PrintInputReport();

        // Create a circulant matrix
        DistMatrix<Complex<double> > A;
        std::vector<Complex<double> > a( n );
        for( int j=0; j<n; ++j )
            a[j] = j;
        Circulant( A, a );
        if( print )
            A.Print("Circulant matrix:");
#ifdef HAVE_QT5
        if( display )
            Display( A, "Circulant" );
#endif

        // Create a Fourier matrix, which can be used to diagonalize circulant
        // matrices
        DistMatrix<Complex<double> > F;
        Fourier( F, n );
        if( print )
            F.Print("DFT matrix:");
#ifdef HAVE_QT5
        if( display )
            Display( F, "DFT matrix" );
#endif
        
        // Form B := A F
        DistMatrix<Complex<double> > B;
        Zeros( B, n, n );
        Gemm( NORMAL, NORMAL, 
              Complex<double>(1), A, F, Complex<double>(0), B );

        // Form A := F^H B = F^H \hat A F
        Gemm( ADJOINT, NORMAL,
              Complex<double>(1), F, B, Complex<double>(0), A );
        if( print )
            A.Print("A := F^H A F");
#ifdef HAVE_QT5
        if( display )
            Display( A, "F^H A F" );
#endif

        // Form the thresholded result
        const int localHeight = A.LocalHeight();
        const int localWidth = A.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                const double absValue = Abs(A.GetLocal(iLocal,jLocal));
                if( absValue < 1e-13 )
                    A.SetLocal(iLocal,jLocal,0);
            }
        }
        if( print )
            A.Print("A with values below 1e-13 removed");
#ifdef HAVE_QT5
        if( display )
            Display( A, "Thresholded (1e-13) A" );
#endif
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
