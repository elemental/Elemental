/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

template<typename T>
void TimedGemm
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C );

int main( int argc, char *argv[] ) 
{
    Initialize( argc, argv );
    const mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );
    const Int commSize = mpi::CommSize( comm );

    try 
    {
        const Int m = Input("--m","height of C",1000);
        const Int n = Input("--n","width of C",1000);     
        const Int k = Input("--k","inner dimension",1000);
        const Int nb = Input("--nb","algorithmic blocksize",128);
        const bool inst = Input("--inst","instrument Gemm?",true);
        Int r = Input("--r","process grid height",0);
        ProcessInput();

        SetBlocksize( nb );
        // If no process grid height was specified, try for a square
        if( r == 0 )
            r = Grid::FindFactor( commSize );
        Grid g( MPI_COMM_WORLD, r );
        if( commRank == 0 )
            std::cout << "g: " << g.Height() << " x " << g.Width() << std::endl;

        Matrix<double> A, B, C;
        Uniform( A, m, k );
        Uniform( B, k, n );
        Uniform( C, m, n );
        mpi::Barrier( comm );

        Timer timer;
        if( commRank == 0 ) 
        {
            timer.Start();
            Gemm( NORMAL, NORMAL, 1., A, B, C );    
            const double gemmTime = timer.Stop();
            const double gflops = (2.*m*n*k)/(gemmTime*1.e9);
            std::cout << "Sequential: " << gemmTime << " secs and " 
                      << gflops << " GFLops" << std::endl;
        }

        DistMatrix<double,CIRC,CIRC> ARoot(g), BRoot(g); 
        if( commRank == 0 ) 
        {
            timer.Start();
            ARoot.CopyFromRoot( A );
            BRoot.CopyFromRoot( B );
            std::cout << "Populate root node: " << timer.Stop() << " secs"
                      << std::endl;
        }
        else
        {
            ARoot.CopyFromNonRoot();
            BRoot.CopyFromNonRoot();
        }

        mpi::Barrier( comm );
        if( commRank == 0 )
            timer.Start();
        DistMatrix<double> ADist( ARoot ), BDist( BRoot ), CDist(g);
        Zeros( CDist, m, n );
        mpi::Barrier( comm );
        if( commRank == 0 ) 
            std::cout << "Spread from root: " << timer.Stop() << " secs"
                      << std::endl;

        if( commRank == 0 )
            timer.Start();
        if( inst )
            TimedGemm( 1., ADist, BDist, 0., CDist );
        else
            Gemm( NORMAL, NORMAL, 1., ADist, BDist, 0., CDist );    
        mpi::Barrier( comm );
        if( commRank == 0 ) 
            std::cout << "Distributed Gemm: " << timer.Stop() << " secs"
                      << std::endl;

        if( commRank == 0 )
            timer.Start();
        DistMatrix<double,CIRC,CIRC> CRoot( CDist );
        mpi::Barrier( comm );
        if( commRank == 0 )
            std::cout << "Gathered to root: " << timer.Stop() << " secs" 
                      << std::endl;
    } catch( std::exception& e ) { ReportException(e); }
    Finalize();
    return 0;
}

template<typename T>
void TimedGemm
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    Timer timerMC, timerMR, timerGemm;

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1,
                               BB,  B2 );

        //--------------------------------------------------------------------//
        timerMC.Start();
        A1_MC_STAR = A1;
        mpi::Barrier( g.Comm() );
        const double timeMC = timerMC.Stop();
        if( g.Rank() == 0 )
        {
            const Int mLocal = A1_MC_STAR.LocalHeight();
            const Int nLocal = A1_MC_STAR.LocalWidth();
            const double mbps = (1.*mLocal*nLocal*sizeof(T))/(timeMC*1.e6);
            std::cout << "[MC,* ] AllGather: " << timeMC
                      << " secs, " << mbps << " MB/s" << " for "
                      << mLocal << " x " << nLocal << " local matrix"
                      << std::endl;
        }
        timerMR.Start();
        B1Trans_MR_STAR.TransposeFrom( B1 );
        mpi::Barrier( g.Comm() );
        const double timeMR = timerMR.Stop();
        if( g.Rank() == 0 )
        {
            const Int nLocal = B1Trans_MR_STAR.LocalHeight();
            const Int mLocal = B1Trans_MR_STAR.LocalWidth();
            const double mbps = (1.*mLocal*nLocal*sizeof(T))/(timeMR*1.e6);
            std::cout << "[* ,MR] AllGather: " << timeMR
                      << " secs, " << mbps << " MB/s" << " for "
                      << mLocal << " x " << nLocal << " local matrix"
                      << std::endl;
        }

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        timerGemm.Start();
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, A1_MC_STAR, B1Trans_MR_STAR, T(1), C );
        mpi::Barrier( g.Comm() );
        const double gemmTime = timerGemm.Stop();
        if( g.Rank() == 0 )
        {
            const Int mLocal = C.LocalHeight();
            const Int nLocal = C.LocalWidth();
            const Int kLocal = A1_MC_STAR.LocalWidth();
            const double gflops = (2.*mLocal*nLocal*kLocal)/(gemmTime*1.e9);
            std::cout << "Local gemm: " << gemmTime << " secs and "
                      << gflops << " GFlops for " << mLocal << " x " << nLocal
                      << " x " << kLocal << " product" << std::endl;
        }
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
}
