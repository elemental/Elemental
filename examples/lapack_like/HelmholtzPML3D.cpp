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

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    typedef double Real;
    typedef Complex<Real> C;

    try
    {
        const Int n1 = Input("--n1","first grid dimension",30);
        const Int n2 = Input("--n2","second grid dimension",30);
        const Int n3 = Input("--n3","third grid dimension",30);
        const double omega = Input("--omega","angular frequency",18.);
        const Int b = Input("--pmlWidth","number of grid points of PML",5);
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
        const Int cutoff = Input("--cutoff","cutoff for nested dissection",128);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();

        const Int N = n1*n2*n3;
        DistSparseMatrix<C> A( N, comm );
        HelmholtzPML( A, n1, n2, n3, C(omega), b, sigma, p );
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        if( commRank == 0 )
            cout << "Generating point-source for y..." << endl;
        DistMultiVec<C> y( N, 1, comm ), z( N, 1, comm );
        Zero( z );
        const Int xSource = n1/2;
        const Int ySource = n2/2;
        const Int zSource = n3/2;
        const Int iSource = xSource + ySource*n1 + zSource*n1*n2;
        const Int firstLocalRow = z.FirstLocalRow();
        const Int localHeight = z.LocalHeight();
        if( iSource >= firstLocalRow && iSource < firstLocalRow+localHeight )
            z.SetLocal( iSource-firstLocalRow, 0, Complex<double>(1.0,0.0) );
        y = z;

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
            cout << "Running block LDL^T..." << endl;
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
                cout << "done, " << mpi::Time()-svdStart << " seconds\n"
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
            mpi::Barrier( comm );
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
            cout << "Checking residual norm of solution..." << endl;
        const double bNorm = Nrm2( z );
        Multiply( NORMAL, C(-1), A, y, C(1), z );
        const double errorNorm = Nrm2( z );
        if( commRank == 0 )
            cout << "|| b     ||_2 = " << bNorm << "\n"
                 << "|| error ||_2 / || b ||_2 = " 
                 << errorNorm/bNorm << "\n" << endl;
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
