/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
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
        ProcessInput();
        PrintInputReport();

        // Create a circulant matrix
        DistMatrix<Complex<double> > A;
        std::vector<Complex<double> > a( n );
        for( int j=0; j<n; ++j )
            a[j] = j;
        Circulant( a, A );
        if( print )
            A.Print("Circulant matrix:");

        // Create a discrete Fourier matrix, which can be used to diagonalize
        // circulant matrices
        DistMatrix<Complex<double> > F;
        DiscreteFourier( n, F );
        if( print )
            F.Print("DFT matrix (F):");
        
        // Form B := A F
        DistMatrix<Complex<double> > B;
        Zeros( n, n, B );
        Gemm( NORMAL, NORMAL, 
              Complex<double>(1), A, F, Complex<double>(0), B );

        // Form A := F^H B = F^H \hat A F
        Gemm( ADJOINT, NORMAL,
              Complex<double>(1), F, B, Complex<double>(0), A );
        if( print )
            A.Print("A := F^H A F");

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
