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

Complex<double> PML( double x, double w, double p, double sigma, double k )
{
    DEBUG_ONLY(
        if( x < 0 || x > w+1e-10 )
            LogicError("Evaluation point not in PML interval");
    )
    const double realPart = 1.0;
    const double arg = x/w;
    const double imagPart = (sigma/w)*std::pow(arg,p)/k;
    return Complex<double>(realPart,imagPart); 
}

Complex<double> 
sInv( int j, int n, int b, double h, double p, double sigma, double k )
{
    if( j < b-1 )
        return PML( (b-1-j)*h, b*h, p, sigma, k );
    else if( j > n-b )
        return PML( (j-(n-b))*h, b*h, p, sigma, k );
    else
        return Complex<double>(1.0,0.0);
}

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
        const double L1 = Input("--L1","length of domain in first dir",1.);
        const double L2 = Input("--L2","length of domain in second dir",1.);
        const double L3 = Input("--L3","length of domain in third dir",1.);
        const int b = Input("--pmlWidth","number of grid points of PML",5);
        const double sigma = Input("--sigma","magnitude of PML profile",1.5);
        const double p = Input("--exponent","exponent of PML profile",3.);
        const bool selInv = Input("--selInv","selectively invert?",false);
        const bool intraPiv = Input("--intraPiv","frontal pivoting?",false);
        const bool natural = Input("--natural","analytical partitions?",true);
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

        const double k = omega/(2*M_PI);
        const int N = n1*n2*n3;
        DistSparseMatrix<C> A( N, comm );
        const double h1 = L1/(n1+1);
        const double h2 = L2/(n2+1);
        const double h3 = L3/(n3+1);
        const double h1Squared = h1*h1;
        const double h2Squared = h2*h2;
        const double h3Squared = h3*h3;

        // Fill our portion of the 3D Helmholtz operator 
        if( commRank == 0 )
        {
            std::cout << "Filling local portion of matrix...";
            std::cout.flush();
        }
        const double fillStart = mpi::Time();
        const int firstLocalRow = A.FirstLocalRow();
        const int localHeight = A.LocalHeight();
        A.StartAssembly();
        A.Reserve( 7*localHeight );
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = firstLocalRow + iLocal;
            const int x = i % n1;
            const int y = (i/n1) % n2;
            const int z = i/(n1*n2);

            const C s1InvL = sInv( x-1, n1, b, h1, p, sigma, k );
            const C s1InvM = sInv( x,   n1, b, h1, p, sigma, k );
            const C s1InvR = sInv( x+1, n1, b, h1, p, sigma, k );
            const C s2InvL = sInv( y-1, n2, b, h2, p, sigma, k );
            const C s2InvM = sInv( y,   n2, b, h2, p, sigma, k );
            const C s2InvR = sInv( y+1, n2, b, h2, p, sigma, k );
            const C s3InvL = sInv( z-1, n3, b, h3, p, sigma, k );
            const C s3InvM = sInv( z,   n3, b, h3, p, sigma, k );
            const C s3InvR = sInv( z+1, n3, b, h3, p, sigma, k );

            const C xTop = s2InvM*s3InvM;
            const C xTempL = xTop/s1InvL;
            const C xTempM = xTop/s1InvM;
            const C xTempR = xTop/s1InvR;
            const C xTermL = (xTempL+xTempM) / (2*h1Squared);
            const C xTermR = (xTempM+xTempR) / (2*h1Squared);

            const C yTop = s1InvM*s3InvM;
            const C yTempL = yTop/s2InvL;
            const C yTempM = yTop/s2InvM;
            const C yTempR = yTop/s2InvR;
            const C yTermL = (yTempL+yTempM) / (2*h2Squared);
            const C yTermR = (yTempM+yTempR) / (2*h2Squared);

            const C zTop = s1InvM*s2InvM;
            const C zTempL = zTop/s3InvL;
            const C zTempM = zTop/s3InvM;
            const C zTempR = zTop/s3InvR;
            const C zTermL = (zTempL+zTempM) / (2*h3Squared);
            const C zTermR = (zTempM+zTempR) / (2*h3Squared);

            const C mainTerm = (xTermL+xTermR+yTermL+yTermR+zTermL+zTermR)
                - omega*omega*s1InvM*s2InvM*s3InvM;

            A.Update( i, i, mainTerm );
            if( x != 0 )
                A.Update( i, i-1, -xTermL );
            if( x != n1-1 )
                A.Update( i, i+1, -xTermR );
            if( y != 0 )
                A.Update( i, i-n1, -yTermL );
            if( y != n2-1 )
                A.Update( i, i+n1, -yTermR );
            if( z != 0 )
                A.Update( i, i-n1*n2, -zTermL );
            if( z != n3-1 )
                A.Update( i, i+n1*n2, -zTermR );
        } 
        A.StopAssembly();
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
            std::cout << "Generating point-source for y..." << std::endl;
        DistMultiVec<C> y( N, 1, comm ), z( N, 1, comm );
        Zero( z );
        const int xSource = n1/2;
        const int ySource = n2/2;
        const int zSource = n3/2;
        const int iSource = xSource + ySource*n1 + zSource*n1*n2;
        if( iSource >= firstLocalRow && iSource < firstLocalRow+localHeight )
            z.SetLocal( iSource-firstLocalRow, 0, Complex<double>(1.0,0.0) );
        y = z;

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
            NaturalNestedDissection
            ( n1, n2, n3, graph, map, sepTree, info, cutoff );
        else
#ifdef HAVE_PARMETIS
            NestedDissection
            ( graph, map, sepTree, info, 
              sequential, numDistSeps, numSeqSeps, cutoff );
#else 
            std::cout << "NestedDissection requested but not built with parmetis" << std::endl;
            return 1;
#endif 
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
            mpi::Barrier( comm );
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
            std::cout << "Checking residual norm of solution..." << std::endl;
        const double bNorm = Norm( z );
        Multiply( C(-1), A, y, C(1), z );
        const double errorNorm = Norm( z );
        if( commRank == 0 )
        {
            std::cout << "|| b     ||_2 = " << bNorm << "\n"
                      << "|| error ||_2 / || b ||_2 = " 
                      << errorNorm/bNorm << "\n"
                      << std::endl;
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
