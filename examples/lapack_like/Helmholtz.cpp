/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int commRank = El::mpi::Rank( comm );

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;
        El::Timer timer;

        const El::Int n1 = El::Input("--n1","first grid dimension",30);
        const El::Int n2 = El::Input("--n2","second grid dimension",30);
        const El::Int n3 = El::Input("--n3","third grid dimension",30);
        const Real omega = El::Input("--omega","angular frequency",Real(18));
        const Real damping = El::Input("--damping","damping parameter",Real(7));
        const bool selInv = El::Input("--selInv","selectively invert?",false);
        const bool intraPiv = El::Input("--intraPiv","frontal pivoting?",false);
        const bool natural = El::Input("--natural","analytic partitions?",true);
        const bool sequential = El::Input
          ("--sequential","sequential partitions?",true);
        const int numDistSeps = El::Input
          ("--numDistSeps",
           "number of separators to try per distributed partition",1);
        const int numSeqSeps = El::Input
          ("--numSeqSeps",
           "number of separators to try per sequential partition",1);
        const int cutoff =
          El::Input("--cutoff","cutoff for nested dissection",128);
        const bool print = El::Input("--print","print matrix?",false);
        const bool display = El::Input("--display","display matrix?",false);
        El::ProcessInput();

        const El::Grid grid(comm);
        El::DistSparseMatrix<Scalar> A(grid);
        Scalar dampedOmega( omega, damping );
        El::Helmholtz( A, n1, n2, n3, dampedOmega*dampedOmega );
        const El::Int N = A.Height();
        if( display )
            El::Display( A, "A" );
        if( print )
            El::Print( A, "A" );

        if( commRank == 0 )
            El::Output("Generating random vector x and forming y := A x");
        timer.Start();
        El::DistMultiVec<Scalar> x( N, 1, grid ), y( N, 1, grid );
        El::MakeUniform( x );
        El::Zero( y );
        El::Multiply( El::NORMAL, Scalar(1), A, x, Scalar(0), y );
        const Real yOrigNorm = El::FrobeniusNorm( y );
        El::mpi::Barrier(comm);
        if( commRank == 0 )
            El::Output(timer.Stop()," seconds");

        El::BisectCtrl ctrl;
        ctrl.sequential = sequential;
        ctrl.numSeqSeps = numSeqSeps;
        ctrl.numDistSeps = numDistSeps;
        ctrl.cutoff = cutoff;

        El::DistSparseLDLFactorization<Scalar> sparseLDLFact;
        const bool hermitian = false;
        if( commRank == 0 )
            El::Output("Running analysis...");
        timer.Start();
        if( natural )
        {
            sparseLDLFact.Initialize3DGridGraph
            ( n1, n2, n3, A, hermitian, ctrl );
        }
        else
        {
            sparseLDLFact.Initialize( A, hermitian, ctrl );
        }
        El::mpi::Barrier(comm);
        if( commRank == 0 )
            El::Output(timer.Stop()," seconds");

        auto& front = sparseLDLFact.Front();
        auto& info = sparseLDLFact.NodeInfo();
        const El::Int rootSepSize = info.size;
        if( commRank == 0 )
            El::Output(rootSepSize," vertices in root separator\n");

        if( commRank == 0 )
            El::Output("Running LDL factorization...");
        El::mpi::Barrier(comm);
        timer.Start();
        El::LDLFrontType type;
        if( intraPiv )
            type = selInv ? El::LDL_INTRAPIV_SELINV_2D : El::LDL_INTRAPIV_2D;
        else
            type = selInv ? El::LDL_SELINV_2D : El::LDL_2D;
        sparseLDLFact.Factor( type );
        El::mpi::Barrier(comm);
        if( commRank == 0 )
            El::Output(timer.Stop()," seconds");

        if( info.child != nullptr && info.child->onLeft )
        {
            if( commRank == 0 )
                El::Output
                ("Computing SVD of connectivity of second separator to "
                 "the root separator...");
            timer.Start();
            const auto& FL = front.child->L2D;
            const El::Grid& grid = FL.Grid();
            const El::Int height = FL.Height();
            const El::Int width = FL.Width();
            auto B = FL( El::IR(width,height), El::IR(0,width) );
            El::DistMatrix<Real,El::VR,El::STAR> singVals_VR_STAR( grid );
            El::SVD( B, singVals_VR_STAR );
            El::DistMatrix<Real,El::CIRC,El::CIRC> singVals( singVals_VR_STAR );
            El::mpi::Barrier( grid.Comm() );
            const Real twoNorm = El::MaxNorm( singVals_VR_STAR );
            const El::Int minDim = singVals_VR_STAR.Height();
            if( grid.Rank() == singVals.Root() )
            {
                El::Output
                (" two-norm is ",twoNorm," (took ",timer.Stop()," seconds)");
                for( Real tol=1e-1; tol>=Real(1e-10); tol/=10 )
                {
                    El::Int numRank = minDim;
                    for( El::Int j=0; j<minDim; ++j )
                    {
                        if( singVals.GetLocal(j,0) <= twoNorm*tol )
                        {
                            numRank = j;
                            break;
                        }
                    }
                    El::Output("  rank (",tol,")=",numRank,"/",minDim);
                }
            }
        }

        if( commRank == 0 )
            El::Output
            ("Computing SVD of the largest off-diagonal block of "
             "numerical Green's function on root separator...");
        {
            timer.Start();
            const auto& FL = front.L2D;
            const El::Grid& grid = FL.Grid();
            const El::Int lHalf = rootSepSize/2;
            const El::Int uHalf = rootSepSize - lHalf;
            if( commRank == 0 )
                El::Output("lower half=",lHalf,", upper half=",uHalf);
            auto offDiagBlock =
              FL( El::IR(lHalf,rootSepSize), El::IR(0,lHalf) );
            El::DistMatrix<Real,El::VR,El::STAR> singVals_VR_STAR( grid );
            El::SVD( offDiagBlock, singVals_VR_STAR );
            El::DistMatrix<Real,El::CIRC,El::CIRC> singVals( singVals_VR_STAR );
            El::mpi::Barrier( grid.Comm() );
            const Real twoNorm = El::MaxNorm( singVals_VR_STAR );
            if( grid.Rank() == singVals.Root() )
            {
                El::Output(timer.Stop()," seconds");
                for( Real tol=1e-1; tol>=Real(1e-10); tol/=10 )
                {
                    El::Int numRank = lHalf;
                    for( El::Int j=0; j<lHalf; ++j )
                    {
                        if( singVals.GetLocal(j,0) <= twoNorm*tol )
                        {
                            numRank = j;
                            break;
                        }
                    }
                    El::Output("  rank (",tol,")=",numRank,"/",lHalf);
                }
            }
        }

        if( commRank == 0 )
            El::Output("Solving against y...");
        timer.Start();
        sparseLDLFact.Solve( y );
        El::mpi::Barrier(comm);
        if( commRank == 0 )
            El::Output(timer.Stop()," seconds");

        if( commRank == 0 )
            El::Output("Checking error in computed solution...");
        const Real xNorm = El::FrobeniusNorm( x );
        const Real yNorm = El::FrobeniusNorm( y );
        y -= x;
        const Real errorNorm = El::FrobeniusNorm( y );
        if( commRank == 0 )
            El::Output
            ("|| x     ||_2 = ",xNorm,"\n",
             "|| xComp ||_2 = ",yNorm,"\n",
             "|| A x   ||_2 = ",yOrigNorm,"\n",
             "|| error ||_2 / || x ||_2 = ",errorNorm/xNorm,"\n",
             "|| error ||_2 / || A x ||_2 = ",errorNorm/yOrigNorm);
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
