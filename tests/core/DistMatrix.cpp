/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename T,Dist AColDist,Dist ARowDist,Dist BColDist,Dist BRowDist>
void
Check( DistMatrix<T,AColDist,ARowDist>& A, 
       DistMatrix<T,BColDist,BRowDist>& B, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("Check"))
    const Grid& g = A.Grid();

    const Int height = B.Height();
    const Int width = B.Width();

    OutputFromRoot
    (g.Comm(),
     "Testing [",DistToString(AColDist),",",DistToString(ARowDist),"]",
     " <- [",DistToString(BColDist),",",DistToString(BRowDist),"]");
    Int colAlign = SampleUniform<Int>(0,A.ColStride());
    Int rowAlign = SampleUniform<Int>(0,A.RowStride());
    mpi::Broadcast( colAlign, 0, g.Comm() );
    mpi::Broadcast( rowAlign, 0, g.Comm() );
    A.Align( colAlign, rowAlign );
    A = B;
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError
        ("A ~ ",A.Height()," x ",A.Width(),", B ~ ",B.Height()," x ",B.Width());

    DistMatrix<T,STAR,STAR> A_STAR_STAR(A), B_STAR_STAR(B);
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
        OutputFromRoot(g.Comm(),"PASSED");
        if( print )
            Print( A, "A" );
        if( print ) 
            Print( B, "B" );
    }
    else
    {
        OutputFromRoot(g.Comm(),"FAILED");
        if( print )
            Print( A, "A" );
        if( print ) 
            Print( B, "B" );
        LogicError("Redistribution test failed");
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

    const T center = 0;
    const Base<T> radius = 5;
    Uniform( A, m, n, center, radius );

    {
      DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC(g);
      Check( A_CIRC_CIRC, A, print );
    }

    {
      DistMatrix<T,MC,MR> A_MC_MR(g);
      Check( A_MC_MR, A, print );
    }

    {
      DistMatrix<T,MC,STAR> A_MC_STAR(g);
      Check( A_MC_STAR, A, print );
    }

    {
      DistMatrix<T,MD,STAR> A_MD_STAR(g);
      Check( A_MD_STAR, A, print );
    }

    {
      DistMatrix<T,MR,MC> A_MR_MC(g);
      Check( A_MR_MC, A, print );
    }

    {
      DistMatrix<T,MR,STAR> A_MR_STAR(g);
      Check( A_MR_STAR, A, print );
    }

    {
      DistMatrix<T,STAR,MC> A_STAR_MC(g);
      Check( A_STAR_MC, A, print );
    }

    {
      DistMatrix<T,STAR,MD> A_STAR_MD(g);
      Check( A_STAR_MD, A, print );
    }

    {
      DistMatrix<T,STAR,MR> A_STAR_MR(g);
      Check( A_STAR_MR, A, print );
    }

    {
      DistMatrix<T,STAR,STAR> A_STAR_STAR(g);
      Check( A_STAR_STAR, A, print );
    }

    {
      DistMatrix<T,STAR,VC> A_STAR_VC(g);
      Check( A_STAR_VC, A, print );
    }

    {
      DistMatrix<T,STAR,VR> A_STAR_VR(g);
      Check( A_STAR_VR, A, print );
    }

    {
      DistMatrix<T,VC,STAR> A_VC_STAR(g);
      Check( A_VC_STAR, A, print );
    }

    {
      DistMatrix<T,VR,STAR> A_VR_STAR(g);
      Check( A_VR_STAR, A, print );
    }
}

template<typename T>
void
DistMatrixTest( Int m, Int n, const Grid& g, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixTest"))
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<T>());
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
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",50);
        const Int n = Input("--width","width of matrix",50);
        const bool print = Input("--print","print wrong matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );

        DistMatrixTest<Int>( m, n, g, print );

        DistMatrixTest<float>( m, n, g, print );
        DistMatrixTest<Complex<float>>( m, n, g, print );

        DistMatrixTest<double>( m, n, g, print );
        DistMatrixTest<Complex<double>>( m, n, g, print );

#ifdef EL_HAVE_QD
        DistMatrixTest<DoubleDouble>( m, n, g, print );
        DistMatrixTest<QuadDouble>( m, n, g, print );
#endif

#ifdef EL_HAVE_QUAD
        DistMatrixTest<Quad>( m, n, g, print );
        DistMatrixTest<Complex<Quad>>( m, n, g, print );
#endif

#ifdef EL_HAVE_MPC
        DistMatrixTest<BigInt>( m, n, g, print );
        OutputFromRoot(comm,"Setting BigInt precision to 512 bits");
        mpfr::SetMinIntBits( 512 );
        DistMatrixTest<BigInt>( m, n, g, print );

        DistMatrixTest<BigFloat>( m, n, g, print );
        DistMatrixTest<Complex<BigFloat>>( m, n, g, print );
        OutputFromRoot(comm,"Setting BigFloat precision to 512 bits");
        mpfr::SetPrecision( 512 );
        DistMatrixTest<BigFloat>( m, n, g, print );
        DistMatrixTest<Complex<BigFloat>>( m, n, g, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
