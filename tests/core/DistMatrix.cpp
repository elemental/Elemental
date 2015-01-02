/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename T,Dist AColDist,Dist ARowDist,Dist BColDist,Dist BRowDist>
void
Check( DistMatrix<T,AColDist,ARowDist>& A, 
       DistMatrix<T,BColDist,BRowDist>& B, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("Check"))
    const Grid& g = A.Grid();

    const Int commRank = g.Rank();
    const Int height = B.Height();
    const Int width = B.Width();
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);
    DistMatrix<T,STAR,STAR> B_STAR_STAR(g);

    if( commRank == 0 )
    {
        std::cout << "Testing [" << DistToString(AColDist) << ","
                                 << DistToString(ARowDist) << "]"
                  << " <- ["     << DistToString(BColDist) << ","
                                 << DistToString(BRowDist) << "]...";
        std::cout.flush();
    }
    Int colAlign = SampleUniform<Int>(0,A.ColStride());
    Int rowAlign = SampleUniform<Int>(0,A.RowStride());
    mpi::Broadcast( colAlign, 0, g.Comm() );
    mpi::Broadcast( rowAlign, 0, g.Comm() );
    A.Align( colAlign, rowAlign );
    A = B;

    A_STAR_STAR = A;
    B_STAR_STAR = B;

    Int myErrorFlag = 0;
    for( Int j=0; j<width; ++j )
    {
        for( Int i=0; i<height; ++i )
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

    Int summedErrorFlag;
    mpi::AllReduce( &myErrorFlag, &summedErrorFlag, 1, mpi::SUM, g.Comm() );

    if( summedErrorFlag == 0 )
    {
        if( commRank == 0 )
            std::cout << "PASSED" << std::endl;
    }
    else
    {
        if( commRank == 0 )
            std::cout << "FAILED" << std::endl;
        if( print )
            Print( A, "A" );
        if( print ) 
            Print( B, "B" );
    }
}

template<typename T,Dist U,Dist V>
void CheckAll( Int m, Int n, const Grid& g, bool print )
{
    DistMatrix<T,U,V> A(g);
    Int colAlign = SampleUniform<Int>(0,A.ColStride());
    Int rowAlign = SampleUniform<Int>(0,A.RowStride());
    mpi::Broadcast( colAlign, 0, g.Comm() );
    mpi::Broadcast( rowAlign, 0, g.Comm() );
    A.Align( colAlign, rowAlign );
    Uniform( A, m, n );

    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC(g);
    DistMatrix<T,MC,  MR  > A_MC_MR(g);
    DistMatrix<T,MC,  STAR> A_MC_STAR(g);
    DistMatrix<T,MD,  STAR> A_MD_STAR(g);
    DistMatrix<T,MR,  MC  > A_MR_MC(g);
    DistMatrix<T,MR,  STAR> A_MR_STAR(g);
    DistMatrix<T,STAR,MC  > A_STAR_MC(g);
    DistMatrix<T,STAR,MD  > A_STAR_MD(g);
    DistMatrix<T,STAR,MR  > A_STAR_MR(g);
    DistMatrix<T,STAR,STAR> A_STAR_STAR(g);
    DistMatrix<T,STAR,VC  > A_STAR_VC(g);
    DistMatrix<T,STAR,VR  > A_STAR_VR(g);
    DistMatrix<T,VC,  STAR> A_VC_STAR(g);
    DistMatrix<T,VR,  STAR> A_VR_STAR(g);

    Check( A_CIRC_CIRC, A, print );
    Check( A_MC_MR,     A, print );
    Check( A_MC_STAR,   A, print );
    Check( A_MD_STAR,   A, print );
    Check( A_MR_MC,     A, print );
    Check( A_MR_STAR,   A, print );
    Check( A_STAR_MC,   A, print );
    Check( A_STAR_MD,   A, print );
    Check( A_STAR_MR,   A, print );
    Check( A_STAR_STAR, A, print );
    Check( A_STAR_VC,   A, print );
    Check( A_STAR_VR,   A, print );
    Check( A_VC_STAR,   A, print );
    Check( A_VR_STAR,   A, print );
}

template<typename T>
void
DistMatrixTest( Int m, Int n, const Grid& g, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixTest"))
    CheckAll<T,CIRC,CIRC>( m, n, g, print );
    CheckAll<T,MC,  MR  >( m, n, g, print );
    CheckAll<T,MC,  STAR>( m, n, g, print );
    CheckAll<T,MD,  STAR>( m, n, g, print );
    CheckAll<T,MR,  MC  >( m, n, g, print );
    CheckAll<T,MR,  STAR>( m, n, g, print );
    CheckAll<T,STAR,MC  >( m, n, g, print );
    CheckAll<T,STAR,MD  >( m, n, g, print );
    CheckAll<T,STAR,MR  >( m, n, g, print );
    CheckAll<T,STAR,STAR>( m, n, g, print );
    CheckAll<T,STAR,VC  >( m, n, g, print );
    CheckAll<T,STAR,VR  >( m, n, g, print );
    CheckAll<T,VC,  STAR>( m, n, g, print );
    CheckAll<T,VR,  STAR>( m, n, g, print );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool print = Input("--print","print wrong matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );

        if( commRank == 0 )
            std::cout << "Testing with floats:" << std::endl;
        DistMatrixTest<float>( m, n, g, print );

        if( commRank == 0 )
            std::cout << "Testing with doubles:" << std::endl;
        DistMatrixTest<double>( m, n, g, print );

        if( commRank == 0 )
            std::cout << "Testing with single-precision complex:" << std::endl;
        DistMatrixTest<Complex<float>>( m, n, g, print );
        
        if( commRank == 0 )
            std::cout << "Testing with double-precision complex:" << std::endl;
        DistMatrixTest<Complex<double>>( m, n, g, print );
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
