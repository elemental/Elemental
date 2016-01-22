/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template< typename Matrix>
class AOp {

    public:
    AOp( const Matrix& A_, Orientation o_ = El::NORMAL) : A( A_), o(o_) {}
    
    template< typename V>
    void operator()( V& v) const {
        V w;
        Gemv(o, 1.0, A, v, w);
        v = w;
    }

    Int Height() const {
        if( o == El::NORMAL){ return A.Height(); }
        return A.Width();
    }

    Int Width() const {
        if( o == El::NORMAL){ return A.Width(); }
        return A.Height();
    }

    const Matrix& A;
    const Orientation o = El::NORMAL;
}; //Matrix Operator


template<typename F>
void TestTSVD( const DistMatrix<F>& A, Int k=3){
    typedef Base<F> Real;
    const auto& g = A.Grid();
    int m = A.Height();
    int n = A.Width();
    DistMatrix<F> Uk( g), Sk( g), Vk( g);
    DistMatrix<F> U( g), V( g);
    DistMatrix<Base<F>, VR, STAR> S( g);
    if( g.Rank() == 0 )
        Output("  Starting TSVD factorization...");
    const double startTime = mpi::Time();
    typedef AOp<DistMatrix<F>> Operator; 
    Operator Aop( A);
    Operator AAdjop( A, El::ADJOINT);
    TSVD( m, n, Aop, AAdjop, k, Uk, Sk, Vk);
    const double runTime = mpi::Time() - startTime;
    Output("TSVD Time = ",runTime);
    mpi::Barrier( g.Comm() );
    if( g.Rank() == 0 )
        Output("  Starting SVD factorization..."); 
    SVD( A, U, S, V);  
    mpi::Barrier( g.Comm() );
    for( Int i = 0; i < k; ++i){
        std::cout <<  i << ": ";
        const Real& A = Sk.Get(i,0);
        const Real& B = S.Get(i,0);
        std::cout << S.Get(i,0) << std::endl;
        //std::cout << Norm(Abs(Uk( ALL, IR(i))) - Abs(U( ALL, IR(i)))) << " ";
        //std::cout << Norm(Abs(Vk( ALL, IR(i))) - Abs(V( ALL, IR(i)))) << " ";
    }
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, order );
        SetBlocksize( nb );
        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test TSQR");

        if( commRank == 0 )
            Output("Testing with doubles:");
        DistMatrix<double> A;
        Laplacian(A, m, n);
        TestTSVD( A, 3);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
