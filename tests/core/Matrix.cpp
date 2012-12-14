/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace elem;

template<typename T> 
void TestMatrix( int m, int n, int ldim )
{
    if( m > ldim || ldim == 0 )
        throw std::logic_error("Leading dimension must be >= m and nonzero");
    std::vector<T> buffer(ldim*n);
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            buffer[i+j*ldim] = i+j*m;

    Matrix<T> A( m, n, &buffer[0], ldim );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            if( A.Get(i,j) != buffer[i+j*ldim] )
                throw std::logic_error
                ("Matrix class was not properly filled with buffer");

    const Matrix<T> B( m, n, (const T*)&buffer[0], ldim );

    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            if( B.Get(i,j) != buffer[i+j*ldim] )
                throw std::logic_error
                ("Matrix class was not properly filled with const buffer");

    const int commRank = mpi::CommRank( mpi::COMM_WORLD );
    if( commRank == 0 )
        std::cout << "passed" << std::endl;
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const int ldim = Input("--ldim","leading dimension",100);
        ProcessInput();
        PrintInputReport();

        if( commRank == 0 )
        {
            std::cout << "Testing with doubles...";
            std::cout.flush();
        }
        TestMatrix<double>( m, n, ldim );

        if( commRank == 0 )
        {
            std::cout << "Testing with double-precision complex...";
            std::cout.flush();
        }
        TestMatrix<Complex<double> >( m, n, ldim );
    }
    catch( ArgException& e ) { }
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
