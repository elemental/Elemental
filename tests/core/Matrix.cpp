/*
   Copyright (c) 2009-2011, Jack Poulson
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
using namespace elemental;

void Usage()
{
    cout << "Basic tests of the elemental::Matrix class.\n"
         << "  Matrix <m> <n> <ldim>\n"
         << "    m: height of matrix\n"
         << "    n: width of matrix\n"
         << "    ldim: leading dimension of matrix" << endl;
}

template<typename T> // represents an a real or complex ring
void TestMatrix( int m, int n, int ldim )
{
    if( m > ldim || ldim == 0 )
        throw logic_error( "Leading dimension must be >= m and nonzero." );
    std::vector<T> buffer(ldim*n);
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            buffer[i+j*ldim] = i+j*m;

    Matrix<T> A( m, n, &buffer[0], ldim );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            if( A.Get(i,j) != buffer[i+j*ldim] )
                throw logic_error
                ( "Matrix class was not properly filled with buffer." );

    const Matrix<T> B( m, n, (const T*)&buffer[0], ldim );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            if( B.Get(i,j) != buffer[i+j*ldim] )
                throw logic_error
                ( "Matrix class was not properly filled with const buffer." );

    int rank = mpi::CommRank( mpi::COMM_WORLD );
    if( rank == 0 )
        cout << "passed" << endl;
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    int rank = mpi::CommRank( comm );

    if( argc < 4 )
    {
        if( rank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try 
    {
        int argNum = 0;
        int m = atoi(argv[++argNum]);
        int n = atoi(argv[++argNum]);
        int ldim = atoi(argv[++argNum]);

        if( rank == 0 )
        {
            cout << "Testing with doubles...";
            cout.flush();
        }
        TestMatrix<double>( m, n, ldim );

#ifndef WITHOUT_COMPLEX
        if( rank == 0 )
        {
            cout << "Testing with double-precision complex...";
            cout.flush();
        }
        TestMatrix<dcomplex>( m, n, ldim );
#endif
    }
    catch( exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        cerr << "Process " << rank << " caught error message:\n" 
             << e.what() << endl;
    }   
    Finalize();
    return 0;
}

