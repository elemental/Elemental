/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

void Usage()
{
    cout << "SVD <m> <n>\n"
         << "  <m>: height of random matrix to test SVD on\n"
         << "  <n>: width of random matrix to test SVD on\n"
         << endl;
}

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 3 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }
    const int m = atoi( argv[1] );
    const int n = atoi( argv[2] );

    try 
    {
        Grid g( comm );
        DistMatrix<C> A( g );
        Uniform( m, n, A );

        // Compute just the singular values 
        DistMatrix<R,VR,STAR> sOnly( g );
        DistMatrix<C> U( A );
        SingularValues( U, sOnly );

        // Compute the SVD of A 
        DistMatrix<C> V( g );
        DistMatrix<R,VR,STAR> s( g );
        U = A;
        SVD( U, s, V );

        // Compare the singular values from both methods
        Axpy( (R)-1, s, sOnly );
        const R singValDiff = Norm( sOnly, FROBENIUS_NORM );

        const R twoNormOfA = Norm( s, MAX_NORM );
        const R maxNormOfA = Norm( A, MAX_NORM );
        const R oneNormOfA = Norm( A, ONE_NORM );
        const R infNormOfA = Norm( A, INFINITY_NORM );
        const R frobNormOfA = Norm( A, FROBENIUS_NORM );

        DiagonalScale( RIGHT, NORMAL, s, U );
        Gemm( NORMAL, ADJOINT, (C)-1, U, V, (C)1, A );
        const R maxNormOfE = Norm( A, MAX_NORM );
        const R oneNormOfE = Norm( A, ONE_NORM );
        const R infNormOfE = Norm( A, INFINITY_NORM );
        const R frobNormOfE = Norm( A, FROBENIUS_NORM );
        const R epsilon = lapack::MachineEpsilon<R>();
        const R scaledResidual = frobNormOfE / (max(m,n)*epsilon*twoNormOfA);

        if( commRank == 0 )
        {
            cout << "||A||_max   = " << maxNormOfA << "\n"
                 << "||A||_1     = " << oneNormOfA << "\n"
                 << "||A||_oo    = " << infNormOfA << "\n"
                 << "||A||_F     = " << frobNormOfA << "\n"
                 << "||A||_2     = " << twoNormOfA << "\n"
                 << "\n"
                 << "||A - U Sigma V^H||_max = " << maxNormOfE << "\n"
                 << "||A - U Sigma V^H||_1   = " << oneNormOfE << "\n"
                 << "||A - U Sigma V^H||_oo  = " << infNormOfE << "\n"
                 << "||A - U Sigma V^H||_F   = " << frobNormOfE << "\n"
                 << "||A - U Sigma V_H||_F / (max(m,n) eps ||A||_2) = " 
                 << scaledResidual << "\n" 
                 << "\n"
                 << "|| sError ||_2 = " << singValDiff << std::endl;
        }
    }
    catch( exception& e )
    {
        cerr << "Process " << commRank << " caught exception with message: "
             << e.what() << endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}

