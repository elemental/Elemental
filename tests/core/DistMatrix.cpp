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
#include <cstdlib>
#include <ctime>
#include "elemental.hpp"
using namespace elem;

void Usage()
{
    std::cout 
        << "Test all of the aligned redistributions of DistMatrix class.\n\n"
        << "  DistMatrix <r> <c> <m> <n>\n\n"
        << "  r: number of process rows\n"
        << "  c: number of process cols\n"
        << "  m: height of matrices\n"
        << "  n: width of matrices\n" << std::endl;
}

// represents a real or complex ring
template<typename T, Distribution AColDist, Distribution ARowDist,
                     Distribution BColDist, Distribution BRowDist>
void
Check( DistMatrix<T,AColDist,ARowDist>& A, 
       DistMatrix<T,BColDist,BRowDist>& B )
{
#ifndef RELEASE
    PushCallStack("Check");
#endif
    const Grid& g = A.Grid();

    const int rank = g.Rank();
    const int height = B.Height();
    const int width = B.Width();
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> B_STAR_STAR(g);

    if( rank == 0 )
    {
        std::cout << "Testing [" << DistToString(AColDist) << ","
                                 << DistToString(ARowDist) << "]"
                  << " <- ["     << DistToString(BColDist) << ","
                                 << DistToString(BRowDist) << "]...";
        std::cout.flush();
    }
    A = B;

    A_STAR_STAR = A;
    B_STAR_STAR = B;

    int myErrorFlag = 0;
    for( int j=0; j<width; ++j )
    {
        for( int i=0; i<height; ++i )
        {
            if( A_STAR_STAR.GetLocal(i,j) != B_STAR_STAR.GetLocal(i,j) )
            {
                myErrorFlag = 1;
                break;
            }
        }
        if( myErrorFlag != 0 )
            break;
    }

    int summedErrorFlag;
    mpi::AllReduce( &myErrorFlag, &summedErrorFlag, 1, mpi::SUM, g.Comm() );

    if( summedErrorFlag == 0 )
        if( rank == 0 )
            std::cout << "PASSED" << std::endl;
    else
        throw std::logic_error("Redistribution failed");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T> // represents a real or complex ring
void
DistMatrixTest( int m, int n, const Grid& g )
{
#ifndef RELEASE
    PushCallStack("DistMatrixTest");
#endif
    DistMatrix<T,MC,  MR  > A_MC_MR(g);
    DistMatrix<T,MC,  STAR> A_MC_STAR(g);
    DistMatrix<T,STAR,MR  > A_STAR_MR(g);
    DistMatrix<T,MR,  MC  > A_MR_MC(g);
    DistMatrix<T,MR,  STAR> A_MR_STAR(g);
    DistMatrix<T,STAR,MC  > A_STAR_MC(g);
    DistMatrix<T,VC,  STAR> A_VC_STAR(g);
    DistMatrix<T,STAR,VC  > A_STAR_VC(g);
    DistMatrix<T,VR,  STAR> A_VR_STAR(g);
    DistMatrix<T,STAR,VR  > A_STAR_VR(g);
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);

    // Communicate from A[MC,MR] 
    Uniform( m, n, A_MC_MR );
    Check( A_MC_STAR,   A_MC_MR );
    Check( A_STAR_MR,   A_MC_MR );
    Check( A_MR_MC,     A_MC_MR );
    Check( A_MR_STAR,   A_MC_MR );
    Check( A_STAR_MC,   A_MC_MR );
    Check( A_VC_STAR,   A_MC_MR );
    Check( A_STAR_VC,   A_MC_MR );
    Check( A_VR_STAR,   A_MC_MR );
    Check( A_STAR_VR,   A_MC_MR );
    Check( A_STAR_STAR, A_MC_MR );

    // Communicate from A[MC,*]
    Uniform( m, n, A_MC_STAR );
    Check( A_MC_MR,     A_MC_STAR );
    Check( A_STAR_MR,   A_MC_STAR );
    Check( A_MR_MC,     A_MC_STAR );
    Check( A_MR_STAR,   A_MC_STAR );
    Check( A_STAR_MC,   A_MC_STAR );
    Check( A_VC_STAR,   A_MC_STAR );
    Check( A_STAR_VC,   A_MC_STAR );
    Check( A_VR_STAR,   A_MC_STAR );
    Check( A_STAR_VR,   A_MC_STAR );
    Check( A_STAR_STAR, A_MC_STAR );

    // Communicate from A[*,MR]
    Uniform( m, n, A_STAR_MR );
    Check( A_MC_MR,     A_STAR_MR );
    Check( A_MC_STAR,   A_STAR_MR );
    Check( A_MR_MC,     A_STAR_MR );
    Check( A_MR_STAR,   A_STAR_MR );
    Check( A_STAR_MC,   A_STAR_MR );
    Check( A_VC_STAR,   A_STAR_MR );
    Check( A_STAR_VC,   A_STAR_MR );
    Check( A_VR_STAR,   A_STAR_MR );
    Check( A_STAR_VR,   A_STAR_MR );
    Check( A_STAR_STAR, A_STAR_MR );
    
    // Communicate from A[MR,MC]
    Uniform( m, n, A_MR_MC );
    Check( A_MC_MR,     A_MR_MC );
    Check( A_MC_STAR,   A_MR_MC );
    Check( A_STAR_MR,   A_MR_MC );
    Check( A_MR_STAR,   A_MR_MC );
    Check( A_STAR_MC,   A_MR_MC );
    Check( A_VC_STAR,   A_MR_MC );
    Check( A_STAR_VC,   A_MR_MC );
    Check( A_VR_STAR,   A_MR_MC );
    Check( A_STAR_VR,   A_MR_MC );
    Check( A_STAR_STAR, A_MR_MC );

    // Communicate from A[MR,*]
    Uniform( m, n, A_MR_STAR );
    Check( A_MC_MR,     A_MR_STAR );
    Check( A_MC_STAR,   A_MR_STAR );
    Check( A_STAR_MR,   A_MR_STAR );
    Check( A_MR_MC,     A_MR_STAR );
    Check( A_STAR_MC,   A_MR_STAR );
    Check( A_VC_STAR,   A_MR_STAR );
    Check( A_STAR_VC,   A_MR_STAR );
    Check( A_VR_STAR,   A_MR_STAR );
    Check( A_STAR_VR,   A_MR_STAR );
    Check( A_STAR_STAR, A_MR_STAR );

    // Communicate from A[*,MC]
    Uniform( m, n, A_STAR_MC );
    Check( A_MC_MR,     A_STAR_MC );
    Check( A_MC_STAR,   A_STAR_MC );
    Check( A_STAR_MR,   A_STAR_MC );
    Check( A_MR_MC,     A_STAR_MC );
    Check( A_MR_STAR,   A_STAR_MC );
    Check( A_VC_STAR,   A_STAR_MC );
    Check( A_STAR_VC,   A_STAR_MC );
    Check( A_VR_STAR,   A_STAR_MC );
    Check( A_STAR_VR,   A_STAR_MC );
    Check( A_STAR_STAR, A_STAR_MC );
 
    // Communicate from A[VC,*]
    Uniform( m, n, A_VC_STAR );
    Check( A_MC_MR,     A_VC_STAR );
    Check( A_MC_STAR,   A_VC_STAR );
    Check( A_STAR_MR,   A_VC_STAR );
    Check( A_MR_MC,     A_VC_STAR );
    Check( A_MR_STAR,   A_VC_STAR );
    Check( A_STAR_MC,   A_VC_STAR );
    Check( A_STAR_VC,   A_VC_STAR );
    Check( A_VR_STAR,   A_VC_STAR );
    Check( A_STAR_VR,   A_VC_STAR );
    Check( A_STAR_STAR, A_VC_STAR );

    // Communicate from A[*,VC]
    Uniform( m, n, A_STAR_VC );
    Check( A_MC_MR,     A_STAR_VC );
    Check( A_MC_STAR,   A_STAR_VC );
    Check( A_STAR_MR,   A_STAR_VC );
    Check( A_MR_MC,     A_STAR_VC );
    Check( A_MR_STAR,   A_STAR_VC );
    Check( A_STAR_MC,   A_STAR_VC );
    Check( A_VC_STAR,   A_STAR_VC );
    Check( A_VR_STAR,   A_STAR_VC );
    Check( A_STAR_VR,   A_STAR_VC );
    Check( A_STAR_STAR, A_STAR_VC );

    // Communicate from A[VR,*]
    Uniform( m, n, A_VR_STAR );
    Check( A_MC_MR,     A_VR_STAR );
    Check( A_MC_STAR,   A_VR_STAR );
    Check( A_STAR_MR,   A_VR_STAR );
    Check( A_MR_MC,     A_VR_STAR );
    Check( A_MR_STAR,   A_VR_STAR );
    Check( A_STAR_MC,   A_VR_STAR );
    Check( A_VC_STAR,   A_VR_STAR );
    Check( A_STAR_VC,   A_VR_STAR );
    Check( A_STAR_VR,   A_VR_STAR );
    Check( A_STAR_STAR, A_VR_STAR );

    // Communicate from A[*,VR]
    Uniform( m, n, A_STAR_VR );
    Check( A_MC_MR,     A_STAR_VR );
    Check( A_MC_STAR,   A_STAR_VR );
    Check( A_STAR_MR,   A_STAR_VR );
    Check( A_MR_MC,     A_STAR_VR );
    Check( A_MR_STAR,   A_STAR_VR );
    Check( A_STAR_MC,   A_STAR_VR );
    Check( A_VC_STAR,   A_STAR_VR );
    Check( A_STAR_VC,   A_STAR_VR );
    Check( A_VR_STAR,   A_STAR_VR );
    Check( A_STAR_STAR, A_STAR_VR );

    // Communicate from A[*,*]
    Uniform( m, n, A_STAR_STAR );
    Check( A_MC_MR,   A_STAR_STAR );
    Check( A_MC_STAR, A_STAR_STAR );
    Check( A_STAR_MR, A_STAR_STAR );
    Check( A_MR_MC,   A_STAR_STAR );
    Check( A_MR_STAR, A_STAR_STAR );
    Check( A_STAR_MC, A_STAR_STAR );
    Check( A_VC_STAR, A_STAR_STAR );
    Check( A_STAR_VC, A_STAR_STAR );
    Check( A_VR_STAR, A_STAR_STAR );
    Check( A_STAR_VR, A_STAR_STAR );
#ifndef RELEASE
    PopCallStack();
#endif
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    if( argc < 5 )
    {
        if( commRank == 0 )
            Usage();
        Finalize();
        return 0;
    }

    try
    {
        int argNum = 0;
        const int r = atoi(argv[++argNum]);
        const int c = atoi(argv[++argNum]);
        const int m = atoi(argv[++argNum]);
        const int n = atoi(argv[++argNum]);
#ifndef RELEASE
        if( commRank == 0 )
        {
            std::cout << "==========================================\n"
                      << " In debug mode! Performance will be poor! \n"
                      << "==========================================" 
                      << std::endl;
        }
#endif
        const Grid g( comm, r, c );

        if( commRank == 0 )
        {
            std::cout << "--------------------\n"
                      << "Testing with floats:\n"
                      << "--------------------" << std::endl;
        }
        DistMatrixTest<float>( m, n, g );

        if( commRank == 0 )
        {
            std::cout << "---------------------\n"
                      << "Testing with doubles:\n"
                      << "---------------------" << std::endl;
        }
        DistMatrixTest<double>( m, n, g );

        if( commRank == 0 )
        {
            std::cout << "--------------------------------------\n"
                      << "Testing with single-precision complex:\n"
                      << "--------------------------------------" << std::endl;
        }
        DistMatrixTest<Complex<float> >( m, n, g );
        
        if( commRank == 0 )
        {
            std::cout << "--------------------------------------\n"
                      << "Testing with double-precision complex:\n"
                      << "--------------------------------------" << std::endl;
        }
        DistMatrixTest<Complex<double> >( m, n, g );
    }
    catch( std::exception& e )
    {
#ifndef RELEASE
        DumpCallStack();
#endif
        std::cerr << "Process " << commRank << " caught error message:\n"
                  << e.what() << std::endl;
    }
    Finalize();
    return 0;
}

