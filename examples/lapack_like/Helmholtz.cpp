/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
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

        DistSparseMatrix<C> A(comm);
        C dampedOmega( omega, damping );
        Helmholtz( A, n1, n2, n3, dampedOmega*dampedOmega );
        const Int N = A.Height();
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        if( commRank == 0 )
            cout << "Generating random vector x and forming y := A x" << endl;
        const double multiplyStart = mpi::Time();
        DistMultiVec<C> x( N, 1, comm ), y( N, 1, comm );
        MakeUniform( x );
        Zero( y );
        Multiply( NORMAL, C(1), A, x, C(0), y );
        const double yOrigNorm = Nrm2( y );
        mpi::Barrier( comm );
        if( commRank == 0 )
            cout << mpi::Time()-multiplyStart << " seconds" << endl;

        if( commRank == 0 )
            cout << "Running nested dissection..." << endl;
        const double nestedStart = mpi::Time();
        const auto& graph = A.DistGraph();
        DistSymmNodeInfo info;
        DistSeparator sep;
        DistMap map, invMap;
        if( natural )
        {
            NaturalNestedDissection
            ( n1, n2, n3, graph, map, sep, info, cutoff );
        }
        else
        { 
            BisectCtrl ctrl;
            ctrl.sequential = sequential;
            ctrl.numSeqSeps = numSeqSeps;
            ctrl.numDistSeps = numDistSeps;
            ctrl.cutoff = cutoff;

            NestedDissection( graph, map, sep, info, ctrl );
        }
        InvertMap( map, invMap );
        mpi::Barrier( comm );
        if( commRank == 0 )
            cout << mpi::Time()-nestedStart << " seconds" << endl;

        const int rootSepSize = info.size;
        if( commRank == 0 )
            cout << rootSepSize << " vertices in root separator\n" << endl;

        if( commRank == 0 )
            cout << "Building DistSymmFront tree..." << endl;
        mpi::Barrier( comm );
        const double buildStart = mpi::Time();
        DistSymmFront<C> front( A, map, sep, info, false );
        mpi::Barrier( comm );
        if( commRank == 0 )
            cout << mpi::Time()-buildStart << " seconds" << endl;

        if( commRank == 0 )
            cout << "Running LDL factorization..." << endl;
        mpi::Barrier( comm );
        const double ldlStart = mpi::Time();
        SymmFrontType type;
        if( intraPiv )
            type = ( selInv ? LDL_INTRAPIV_SELINV_2D : LDL_INTRAPIV_2D );
        else
            type = ( selInv ? LDL_SELINV_2D : LDL_2D );
        LDL( info, front, type );
        mpi::Barrier( comm );
        if( commRank == 0 )
            cout << mpi::Time()-ldlStart << " seconds" << endl;


        if( info.child != nullptr && info.child->onLeft )
        {
            if( commRank == 0 )
                cout << "Computing SVD of connectivity of second separator to "
                        "the root separator..." << endl;
            const double svdStart = mpi::Time();
            const auto& FL = front.child->L2D;
            const Grid& grid = FL.Grid();
            const int height = FL.Height();
            const int width = FL.Width();
            auto B = FL( IR(width,height), IR(0,width) );
            auto BCopy( B );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( BCopy, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            const Int minDim = singVals_VR_STAR.Height();
            if( grid.Rank() == singVals.Root() )
            {
                cout << mpi::Time()-svdStart << " seconds\n"
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
                    cout << "  rank (" << tol << ")=" << numRank 
                         << "/" << minDim << endl;
                }
            }
        }

        if( commRank == 0 )
            cout << "Computing SVD of the largest off-diagonal block of "
                    "numerical Green's function on root separator..." << endl;
        {
            const double svdStart = mpi::Time();
            const auto& FL = front.L2D;
            const Grid& grid = FL.Grid();
            const int lHalf = rootSepSize/2;
            const int uHalf = rootSepSize - lHalf;
            if( commRank == 0 )
                cout << "lower half=" << lHalf
                     << ", upper half=" << uHalf << endl;
            auto offDiagBlock = FL( IR(lHalf,rootSepSize), IR(0,lHalf) );
            auto offDiagBlockCopy( offDiagBlock );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( offDiagBlockCopy, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            if( grid.Rank() == singVals.Root() )
            {
                cout << "done, " << mpi::Time()-svdStart << " seconds\n";
                for( double tol=1e-1; tol>=1e-10; tol/=10 )
                {
                    int numRank = lHalf;
                    for( int j=0; j<lHalf; ++j )
                    {
                        if( singVals.GetLocal(j,0) <= twoNorm*tol )
                        {
                            numRank = j;
                            break;
                        }
                    }
                    cout << "  rank (" << tol << ")=" << numRank
                         << "/" << lHalf << endl;
                }
            }
        }

        if( commRank == 0 )
            cout << "Solving against y..." << endl;
        const double solveStart = mpi::Time();
        ldl::SolveAfter( invMap, info, front, y );
        mpi::Barrier( comm );
        if( commRank == 0 )
            cout << mpi::Time()-solveStart << " seconds" << endl;

        if( commRank == 0 )
            cout << "Checking error in computed solution..." << endl;
        const double xNorm = Nrm2( x );
        const double yNorm = Nrm2( y );
        Axpy( C(-1), x, y );
        const double errorNorm = Nrm2( y );
        if( commRank == 0 )
            cout << "|| x     ||_2 = " << xNorm << "\n"
                 << "|| xComp ||_2 = " << yNorm << "\n"
                 << "|| A x   ||_2 = " << yOrigNorm << "\n"
                 << "|| error ||_2 / || x ||_2 = " 
                 << errorNorm/xNorm << "\n"
                 << "|| error ||_2 / || A x ||_2 = " 
                 << errorNorm/yOrigNorm << endl;
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
