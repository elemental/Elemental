/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/convex/SingularValueSoftThreshold.hpp"
#include "elemental/matrices/Uniform.hpp"
#include <set>
using namespace elem;

// Corrupt a percentage of the entries with uniform samples from the unit ball
template<typename F>
int Corrupt( DistMatrix<F>& A, double percentCorrupt )
{
#ifndef RELEASE
    PushCallStack("Corrupt");
#endif
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int localSize = localHeight*localWidth;
    int numLocalCorrupt = (percentCorrupt/100.)*localSize;
    int numLocalCollisions = 0;
    std::set<int> localIndices;
    for( int k=0; k<numLocalCorrupt; ++k )
    {
        const int localIndex = rand() % localSize;
        if( localIndices.count(localIndex) != 0 )
        {
            ++numLocalCollisions;
            continue;
        } 
        const int iLocal = localIndex % localHeight;
        const int jLocal = localIndex / localHeight;
        const F perturb = SampleUnitBall<F>();
        A.SetLocal( iLocal, jLocal, A.GetLocal(iLocal,jLocal)+perturb );
    }
    
    int numCorrupt;
    numLocalCorrupt -= numLocalCollisions;
    mpi::AllReduce
    ( &numLocalCorrupt, &numCorrupt, 1, mpi::SUM, A.Grid().VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return numCorrupt;
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--height","height of matrix",100);
        const int rank = Input("--rank","rank of structured matrix",10);
        const double percentCorrupt = 
            Input("--percentCorrupt","percentage of corrupted entries",10.);
        const int maxIts = Input("--maxIts","maximum iterations",200);
        const double tau = Input("--tau","step size",1.);
        const double tol = Input("--tol","tolerance",1.e-6);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> LTrue;
        {
            DistMatrix<double> U, V;
            Uniform( n, rank, U );
            Uniform( n, rank, V );
            Zeros( n, n, LTrue );
            const double numCorruptApprox = percentCorrupt*(n/10.)*(n/10.);
            Gemm( NORMAL, ADJOINT, 1./numCorruptApprox, U, V, 0., LTrue );
        }
        const double frobLTrue = FrobeniusNorm( LTrue );
        if( commRank == 0 )
            std::cout << "|| L ||_F = " << frobLTrue << std::endl;
        if( print )
            LTrue.Print("True L");

        DistMatrix<double> STrue;
        Zeros( n, n, STrue );
        const int numCorrupt = Corrupt( STrue, percentCorrupt );
        const double frobNormSTrue = FrobeniusNorm( STrue );
        if( commRank == 0 )
            std::cout << "number of corrupted entries: " << numCorrupt << "\n"
                      << "|| S ||_F = " << frobNormSTrue << std::endl;
        if( print )
            STrue.Print("True S");

        // M = LTrue + STrue
        DistMatrix<double> M( LTrue );
        Axpy( 1., STrue, M );
        const double frobM = FrobeniusNorm( M );
        if( commRank == 0 )
            std::cout << "|| M ||_F = " << frobM << std::endl;

        // Run simple RobustPCA here...
        const double gamma = 1. / sqrt(numCorrupt);
        DistMatrix<double> L, S, Y, E;
        Zeros( n, n, L );
        Zeros( n, n, S );
        Zeros( n, n, Y );
        int numIts = 0;
        while( true )
        {
            // SVT_t(M-S+Y)
            L = M;
            Axpy( -1., S, L );
            Axpy(  1., Y, L );
            // TODO: Modify this routine to return the rank
            SingularValueSoftThreshold( L, tau );
          
            // ST_{gamma*t}(M-L+Y)
            S = M;
            Axpy( -1., L, S );
            Axpy(  1., Y, S );
            SoftThreshold( S, gamma*tau );

            // E := M - (L + S)
            E = M;    
            Axpy( -1., L, E );
            Axpy( -1., S, E );

            const double frobE = FrobeniusNorm( E );
            if( frobE/frobM <= tol )            
            {
                const int numNonzeros = ZeroNorm( S );
                if( commRank == 0 )
                    std::cout << "Converged after " << numIts << " iterations "
                              << " with numNonzeros=" << numNonzeros << " and "
                              << "|| E ||_F / || M ||_F = " << frobE/frobM
                              << std::endl;
                if( print )
                {
                    L.Print("L");
                    S.Print("S"); 
                    E.Print("E");
                }
                break;
            }
            else if( numIts >= maxIts )
            {
                if( commRank == 0 )
                    std::cout << "Aborting after " << maxIts << " iterations"
                              << std::endl;
                break;
            }
            else
            {
                if( commRank == 0 )
                    std::cout << numIts << ": || E ||_F = " << frobE 
                              << std::endl;
            }
            
            // Y := Y + tau E
            Axpy( tau, E, Y );
            ++numIts;
        }
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
