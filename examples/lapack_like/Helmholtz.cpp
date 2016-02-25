/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm;
    const int commRank = mpi::Rank( comm );

    try
    {
        const int n1 = Input("--n1","first grid dimension",30);
        const int n2 = Input("--n2","second grid dimension",30);
        const int n3 = Input("--n3","third grid dimension",30);
        const Real omega = Input("--omega","angular frequency",Real(18));
        const Real damping = Input("--damping","damping parameter",Real(7));
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
            Output("Generating random vector x and forming y := A x");
        const double multiplyStart = mpi::Time();
        DistMultiVec<C> x( N, 1, comm ), y( N, 1, comm );
        MakeUniform( x );
        Zero( y );
        Multiply( NORMAL, C(1), A, x, C(0), y );
        const Real yOrigNorm = Nrm2( y );
        mpi::Barrier();
        if( commRank == 0 )
            Output(mpi::Time()-multiplyStart," seconds");

        if( commRank == 0 )
            Output("Running nested dissection...");
        const double nestedStart = mpi::Time();
        const auto& graph = A.DistGraph();
        ldl::DistNodeInfo info;
        ldl::DistSeparator sep;
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
        mpi::Barrier();
        if( commRank == 0 )
            Output(mpi::Time()-nestedStart," seconds");

        const int rootSepSize = info.size;
        if( commRank == 0 )
            Output(rootSepSize," vertices in root separator\n");

        if( commRank == 0 )
            Output("Building ldl::DistFront tree...");
        mpi::Barrier();
        const double buildStart = mpi::Time();
        ldl::DistFront<C> front( A, map, sep, info, false );
        mpi::Barrier();
        if( commRank == 0 )
            Output(mpi::Time()-buildStart," seconds");

        if( commRank == 0 )
            Output("Running LDL factorization...");
        mpi::Barrier();
        const double ldlStart = mpi::Time();
        LDLFrontType type;
        if( intraPiv )
            type = ( selInv ? LDL_INTRAPIV_SELINV_2D : LDL_INTRAPIV_2D );
        else
            type = ( selInv ? LDL_SELINV_2D : LDL_2D );
        LDL( info, front, type );
        mpi::Barrier();
        if( commRank == 0 )
            Output(mpi::Time()-ldlStart," seconds");


        if( info.child != nullptr && info.child->onLeft )
        {
            if( commRank == 0 )
                Output
                ("Computing SVD of connectivity of second separator to "
                 "the root separator...");
            const double svdStart = mpi::Time();
            const auto& FL = front.child->L2D;
            const Grid& grid = FL.Grid();
            const int height = FL.Height();
            const int width = FL.Width();
            auto B = FL( IR(width,height), IR(0,width) );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( B, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            const Int minDim = singVals_VR_STAR.Height();
            if( grid.Rank() == singVals.Root() )
            {
                Output
                (" two-norm is ",twoNorm," (took ",mpi::Time()-svdStart,
                 " seconds)");
                for( Real tol=1e-1; tol>=Real(1e-10); tol/=10 )
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
                    Output("  rank (",tol,")=",numRank,"/",minDim);
                }
            }
        }

        if( commRank == 0 )
            Output
            ("Computing SVD of the largest off-diagonal block of "
             "numerical Green's function on root separator...");
        {
            const double svdStart = mpi::Time();
            const auto& FL = front.L2D;
            const Grid& grid = FL.Grid();
            const int lHalf = rootSepSize/2;
            const int uHalf = rootSepSize - lHalf;
            if( commRank == 0 )
                Output("lower half=",lHalf,", upper half=",uHalf);
            auto offDiagBlock = FL( IR(lHalf,rootSepSize), IR(0,lHalf) );
            DistMatrix<Real,VR,STAR> singVals_VR_STAR( grid );
            SVD( offDiagBlock, singVals_VR_STAR );
            DistMatrix<Real,CIRC,CIRC> singVals( singVals_VR_STAR );
            mpi::Barrier( grid.Comm() );
            const Real twoNorm = MaxNorm( singVals_VR_STAR );
            if( grid.Rank() == singVals.Root() )
            {
                Output(mpi::Time()-svdStart," seconds");
                for( Real tol=1e-1; tol>=Real(1e-10); tol/=10 )
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
                    Output("  rank (",tol,")=",numRank,"/",lHalf);
                }
            }
        }

        if( commRank == 0 )
            Output("Solving against y...");
        const double solveStart = mpi::Time();
        ldl::SolveAfter( invMap, info, front, y );
        mpi::Barrier();
        if( commRank == 0 )
            Output(mpi::Time()-solveStart," seconds");

        if( commRank == 0 )
            Output("Checking error in computed solution...");
        const Real xNorm = Nrm2( x );
        const Real yNorm = Nrm2( y );
        y -= x;
        const Real errorNorm = Nrm2( y );
        if( commRank == 0 )
            Output
            ("|| x     ||_2 = ",xNorm,"\n",
             "|| xComp ||_2 = ",yNorm,"\n",
             "|| A x   ||_2 = ",yOrigNorm,"\n",
             "|| error ||_2 / || x ||_2 = ",errorNorm/xNorm,"\n",
             "|| error ||_2 / || A x ||_2 = ",errorNorm/yOrigNorm);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
