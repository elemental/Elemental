/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    typedef double Real;
    typedef Complex<Real> C;

    try
    {
        const int n1 = Input("--n1","first grid dimension",30);
        const int n2 = Input("--n2","second grid dimension",30);
        const int n3 = Input("--n3","third grid dimension",30);
        const double omega = Input("--omega","angular frequency",18.);
        const double damping = Input("--damping","damping parameter",7.);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool intraPiv = Input("--intraPiv","frontal pivoting?",false);
        const bool natural = Input("--natural","analytic partitions?",true);
        const bool sequential = Input
            ("--sequential","sequential partitions?",true);
        const int numDistSeps = Input
            ("--numDistSeps",
             "number of separators to try per distributed partition",1);
        const int numSeqSeps = Input
            ("--numSeqSeps",
             "number of separators to try per sequential partition",1);
        const int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        const int N = n1*n2*n3;
        DistSparseMatrix<C> A( N, comm );
        C dampedOmega( omega, damping );
        const double hxInv = n1+1;
        const double hyInv = n2+1;
        const double hzInv = n3+1;
        const double hxInvSquared = hxInv*hxInv;
        const double hyInvSquared = hyInv*hyInv;
        const double hzInvSquared = hzInv*hzInv;
        const C mainTerm = 
            2*(hxInvSquared+hyInvSquared+hzInvSquared) - 
            dampedOmega*dampedOmega;

        // Fill our portion of the 3D Helmholtz operator over the unit-square 
        // using a n1 x n2 x n3 7-point stencil in natural ordering: 
        // (x,y,z) at x + y*n1 + z*n1*n2
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        const double fillStart = mpi::Time();
        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();
        A.Reserve( 7*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n1;
            const int y = (i/n1) % n2;
            const int z = i/(n1*n2);

            A.QueueLocalUpdate( iLocal, i, mainTerm );
            if( x != 0 )
                A.QueueLocalUpdate( iLocal, i-1, -hxInvSquared );
            if( x != n1-1 )
                A.QueueLocalUpdate( iLocal, i+1, -hxInvSquared );
            if( y != 0 )
                A.QueueLocalUpdate( iLocal, i-n1, -hyInvSquared );
            if( y != n2-1 )
                A.QueueLocalUpdate( iLocal, i+n1, -hyInvSquared );
            if( z != 0 )
                A.QueueLocalUpdate( iLocal, i-n1*n2, -hzInvSquared );
            if( z != n3-1 )
                A.QueueLocalUpdate( iLocal, i+n1*n2, -hzInvSquared );
        } 
        A.MakeConsistent();
        mpi::Barrier( comm );
        const double fillStop =  mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << fillStop-fillStart << " seconds" 
                      << std::endl;
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        if( commRank == 0 )
        {
            std::cout << "Generating random vector x and forming y := A x...";
            std::cout.flush();
        }
        const double multiplyStart = mpi::Time();
        DistMultiVec<C> x( N, 1, comm ), y( N, 1, comm );
        MakeUniform( x );
        Zero( y );
        Multiply( NORMAL, C(1), A, x, C(0), y );
        const double yOrigNorm = Nrm2( y );
        mpi::Barrier( comm );
        const double multiplyStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << multiplyStop-multiplyStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running nested dissection...";
            std::cout.flush();
        }
        const double nestedStart = mpi::Time();
        const DistGraph& graph = A.DistGraph();
        DistSymmInfo info;
        DistSeparatorTree sepTree;
        DistMap map, inverseMap;
        if( natural )
        {
            NaturalNestedDissection
            ( n1, n2, n3, graph, map, sepTree, info, cutoff );
        }
        else
        { 
            BisectCtrl ctrl;
            ctrl.sequential = sequential;
            ctrl.numSeqSeps = numSeqSeps;
            ctrl.numDistSeps = numDistSeps;
            ctrl.cutoff = cutoff;

            NestedDissection( graph, map, sepTree, info, ctrl );
        }
        map.FormInverse( inverseMap );
        mpi::Barrier( comm );
        const double nestedStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << nestedStop-nestedStart << " seconds"
                      << std::endl;

        const int rootSepSize = info.distNodes.back().size;
        if( commRank == 0 )
        {
            const int numDistNodes = info.distNodes.size();
            const int numLocalNodes = info.localNodes.size();
            std::cout << "\n"
                      << "On the root process:\n"
                      << "-----------------------------------------\n"
                      << numLocalNodes << " local nodes\n"
                      << numDistNodes  << " distributed nodes\n"
                      << rootSepSize << " vertices in root separator\n"
                      << std::endl;
        }

        if( commRank == 0 )
        {
            std::cout << "Building DistSymmFrontTree...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFrontTree<C> frontTree( A, map, sepTree, info, false );
        mpi::Barrier( comm );
        const double buildStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << buildStop-buildStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Running block LDL^T...";
            std::cout.flush();
        }
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        SymmFrontType frontType;
        if( intraPiv )
            frontType = ( selInv ? LDL_INTRAPIV_SELINV_2D
                                 : LDL_INTRAPIV_2D );
        else
            frontType = ( selInv ? LDL_SELINV_2D
                                 : LDL_2D );
        LDL( info, frontTree, frontType );
        mpi::Barrier( comm );
        const double ldlStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << ldlStop-ldlStart << " seconds" 
                      << std::endl;

        if( commRank == 0 )
        {
            std::cout << "Computing SVD of connectivity of second separator to "
                         "the root separator...";
            std::cout.flush();
        }
        const int numDistFronts = frontTree.distFronts.size();
        if( numDistFronts >= 2 && info.distNodes[numDistFronts-2].onLeft )
        {
            const double svdStart = mpi::Time();
            const DistMatrix<C>& frontL = 
                frontTree.distFronts[numDistFronts-2].front2dL;
            const Grid& grid = frontL.Grid();
            const int height = frontL.Height();
            const int width = frontL.Width();
            auto B = LockedView( frontL, width, 0, height-width, width );
            auto BCopy( B );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( BCopy, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            const Int minDim = singVals_VR_STAR.Height();
            if( grid.Rank() == singVals.Root() )
            {
                std::cout << "done, " << mpi::Time()-svdStart << " seconds\n"
                          << "  two norm=" << twoNorm << "\n";
                for( double tol=1e-1; tol>=1e-10; tol/=10 )
                {
                    int numRank = minDim;
                    for( int j=0; j<minDim; ++j )
                    {
                        if( singVals.GetLocal(j,0) <= twoNorm*tol )
                        {
                            numRank = j;
                            break;
                        }
                    }
                    std::cout << "  rank (" << tol << ")=" << numRank 
                              << "/" << minDim << std::endl;
                }
            }
        }

        if( commRank == 0 )
        {
            std::cout << "Computing SVD of the largest off-diagonal block of "
                         "numerical Green's function on root separator...";
            std::cout.flush();
        }
        {
            const double svdStart = mpi::Time();
            const DistMatrix<C>& front = frontTree.distFronts.back().front2dL;
            const Grid& grid = front.Grid();
            const int lowerHalf = rootSepSize/2;
            const int upperHalf = rootSepSize - lowerHalf;
            if( commRank == 0 )
                std::cout << "lowerHalf=" << lowerHalf
                          << ", upperHalf=" << upperHalf << std::endl;
            auto offDiagBlock = 
                LockedView( front, lowerHalf, 0, upperHalf, lowerHalf );
            auto offDiagBlockCopy( offDiagBlock );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( offDiagBlockCopy, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            if( grid.Rank() == singVals.Root() )
            {
                std::cout << "done, " << mpi::Time()-svdStart << " seconds\n";
                for( double tol=1e-1; tol>=1e-10; tol/=10 )
                {
                    int numRank = lowerHalf;
                    for( int j=0; j<lowerHalf; ++j )
                    {
                        if( singVals.GetLocal(j,0) <= twoNorm*tol )
                        {
                            numRank = j;
                            break;
                        }
                    }
                    std::cout << "  rank (" << tol << ")=" << numRank
                              << "/" << lowerHalf << std::endl;
                }
            }
        }

        if( commRank == 0 )
        {
            std::cout << "Solving against y...";
            std::cout.flush();
        }
        const double solveStart = mpi::Time();
        DistNodalMatrix<C> yNodal;
        yNodal.Pull( inverseMap, info, y );
        Solve( info, frontTree, yNodal );
        yNodal.Push( inverseMap, info, y );
        mpi::Barrier( comm );
        const double solveStop = mpi::Time();
        if( commRank == 0 )
            std::cout << "done, " << solveStop-solveStart << " seconds"
                      << std::endl;

        if( commRank == 0 )
            std::cout << "Checking error in computed solution..." << std::endl;
        const double xNorm = Nrm2( x );
        const double yNorm = Nrm2( y );
        Axpy( C(-1), x, y );
        const double errorNorm = Nrm2( y );
        if( commRank == 0 )
        {
            std::cout << "|| x     ||_2 = " << xNorm << "\n"
                      << "|| xComp ||_2 = " << yNorm << "\n"
                      << "|| A x   ||_2 = " << yOrigNorm << "\n"
                      << "|| error ||_2 / || x ||_2 = " 
                      << errorNorm/xNorm << "\n"
                      << "|| error ||_2 / || A x ||_2 = " 
                      << errorNorm/yOrigNorm
                      << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
