/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename T>
void BasicInstrumentedGemm
( T alpha, const El::DistMatrix<T>& A,
           const El::DistMatrix<T>& B,
  T beta,        El::DistMatrix<T>& C )
{
    const El::Grid& grid = A.Grid();
    C *= beta;

    // Temporary distributions
    El::DistMatrix<T,El::MC,El::STAR> A1_MC_STAR(grid);
    El::DistMatrix<T,El::MR,El::STAR> B1Trans_MR_STAR(grid);
    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    // Start accumulating low-rank updates onto C
    El::Timer timerMC, timerMR, timerGemm;
    const El::Int k = A.Width();
    const El::Int blocksize = El::Blocksize();
    for( El::Int j=0; j<k; j+=blocksize )
    {
        const El::Int nb = El::Min( blocksize, k-j );
        auto A1 = A( El::ALL, El::IR(j,j+nb) );
        auto B1 = B( El::IR(j,j+nb), El::ALL );

        timerMC.Start();
        A1_MC_STAR = A1;
        El::mpi::Barrier( grid.Comm() );
        const double timeMC = timerMC.Stop();
        if( grid.Rank() == 0 )
        {
            const El::Int mLocal = A1_MC_STAR.LocalHeight();
            const El::Int nLocal = A1_MC_STAR.LocalWidth();
            const double mbps = (1.*mLocal*nLocal*sizeof(T))/(timeMC*1.e6);
            El::Output
            ("[MC,* ] AllGather: ",timeMC," secs (",mbps," MB/s) for ",
             mLocal," x ",nLocal," local matrix");
        }
        timerMR.Start();
        El::Transpose( B1, B1Trans_MR_STAR );
        El::mpi::Barrier( grid.Comm() );
        const double timeMR = timerMR.Stop();
        if( grid.Rank() == 0 )
        {
            const El::Int nLocal = B1Trans_MR_STAR.LocalHeight();
            const El::Int mLocal = B1Trans_MR_STAR.LocalWidth();
            const double mbps = (1.*mLocal*nLocal*sizeof(T))/(timeMR*1.e6);
            El::Output
            ("[* ,MR] AllGather: ",timeMR," secs (",mbps," MB/s) for ",
             mLocal," x ",nLocal," local matrix");
        }

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        timerGemm.Start();
        El::LocalGemm
        ( El::NORMAL, El::TRANSPOSE,
          alpha, A1_MC_STAR, B1Trans_MR_STAR, T(1), C );
        El::mpi::Barrier( grid.Comm() );
        const double gemmTime = timerGemm.Stop();
        if( grid.Rank() == 0 )
        {
            const El::Int mLocal = C.LocalHeight();
            const El::Int nLocal = C.LocalWidth();
            const El::Int kLocal = A1_MC_STAR.LocalWidth();
            const double gFlops = (2.*mLocal*nLocal*kLocal)/(gemmTime*1.e9);
            El::Output
            ("Local gemm: ",gemmTime," secs (",gFlops," GFlop/s) for ",mLocal,
             " x ",nLocal," x ",kLocal," product");
        }
    }
}

template<typename T>
void TestGemm
( El::Int m, El::Int n, El::Int k, const El::Grid& grid,
  bool testSequential, bool instrument )
{
    El::Timer timer;

    // Choose arbitrary coefficients.
    const T alpha=2, beta=3;

    if( testSequential && grid.Rank() == 0 )
    {
        El::Matrix<T> A, B, C;
        El::Uniform( A, m, k );
        El::Uniform( B, k, n );
        El::Uniform( C, m, n );

        timer.Start();
        El::Gemm( El::NORMAL, El::NORMAL, alpha, A, B, beta, C );
        const double gemmTime = timer.Stop();
        double gFlops = (2.*m*n*k)/(gemmTime*1.e9);
        if( El::IsComplex<T>::value )
            gFlops *= 4;
        El::Output("Sequential: ",gemmTime," secs (",gFlops," GFlop/s)");
        timer.Start();
    }
    El::mpi::Barrier(grid.Comm());
    if( grid.Rank() == 0 )
    {
        const double rootWaitTime = timer.Stop();
        El::Output("Root waited for ",rootWaitTime," seconds");
    }

    El::DistMatrix<T> A(grid), B(grid), C(grid);
    El::Uniform( A, m, k );
    El::Uniform( B, k, n );
    El::Uniform( C, m, n );

    El::mpi::Barrier(grid.Comm());
    if( grid.Rank() == 0 )
        timer.Start();
    if( instrument )
        BasicInstrumentedGemm( alpha, A, B, beta, C );
    else
        El::Gemm( El::NORMAL, El::NORMAL, alpha, A, B, beta, C );
    El::mpi::Barrier(grid.Comm());
    if( grid.Rank() == 0 )
        El::Output("Distributed Gemm: ",timer.Stop()," secs");
}

int main( int argc, char *argv[] )
{
    El::Environment env( argc, argv );
    const El::mpi::Comm comm;
    const int commRank = El::mpi::Rank( comm );
    const int commSize = El::mpi::Size( comm );

    try
    {
        const El::Int m = El::Input("--m","height of C",1000);
        const El::Int n = El::Input("--n","width of C",1000);
        const El::Int k = El::Input("--k","inner dimension",1000);
        const El::Int blocksize =
          El::Input("--blocksize","algorithmic blocksize",128);
        const bool testSequential =
          El::Input("--testSequential","test sequential Gemm?",true);
        const bool testHigherPrec =
          El::Input("--testHigherPrec","test higher precisions?",false);
        const bool instrument =
          El::Input("--instrument","instrument Gemm?",true);
        El::Int gridHeight = El::Input("--gridHeight","process grid height",0);
        El::ProcessInput();

        El::SetBlocksize( blocksize );

        // If no process grid height was specified, try for a square
        if( gridHeight == 0 )
            gridHeight = El::Grid::DefaultHeight( commSize );
        El::Grid grid( El::mpi::COMM_WORLD, gridHeight );
        if( commRank == 0 )
            El::Output("grid is ",grid.Height()," x ",grid.Width());

        TestGemm<float>
        ( m, n, k, grid, testSequential, instrument );
        TestGemm<El::Complex<float>>
        ( m, n, k, grid, testSequential, instrument );
        TestGemm<double>
        ( m, n, k, grid, testSequential, instrument );
        TestGemm<El::Complex<double>>
        ( m, n, k, grid, testSequential, instrument );

        if( testHigherPrec )
        {
#ifdef EL_HAVE_QD
            TestGemm<El::DoubleDouble>
            ( m, n, k, grid, testSequential, instrument );
            TestGemm<El::Complex<El::DoubleDouble>>
            ( m, n, k, grid, testSequential, instrument );
            TestGemm<El::QuadDouble>
            ( m, n, k, grid, testSequential, instrument );
            TestGemm<El::Complex<El::QuadDouble>>
            ( m, n, k, grid, testSequential, instrument );
#endif
#ifdef EL_HAVE_QUAD
            TestGemm<El::Quad>
            ( m, n, k, grid, testSequential, instrument );
            TestGemm<El::Complex<El::Quad>>
            ( m, n, k, grid, testSequential, instrument );
#endif
#ifdef EL_HAVE_MPC
            TestGemm<El::BigFloat>
            ( m, n, k, grid, testSequential, instrument );
            TestGemm<El::Complex<El::BigFloat>>
            ( m, n, k, grid, testSequential, instrument );
#endif
        }
    } catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
