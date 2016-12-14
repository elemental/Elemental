/*
   Copyright (c) 2009-2016, Jack Poulson
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

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","height of matrix",20);
        const El::Int n = El::Input("--width","width of matrix",100);
        El::Int targetRank = El::Input("--rank","rank of matrix",5);
        El::Int maxSteps = El::Input("--maxSteps","max # of steps of QR",10);
        const Real tol = El::Input("--tol","tolerance for ID",Real(-1));
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );

        const El::Int minDim = El::Min( m, n );
        targetRank = El::Min( targetRank, minDim );
        maxSteps = El::Min( maxSteps, targetRank );

        El::DistMatrix<Scalar> U(grid), V(grid);
        El::Uniform( U, m, targetRank );
        El::Uniform( V, n, targetRank );
        El::DistMatrix<Scalar> A(grid);
        El::Gemm( El::NORMAL, El::ADJOINT, Scalar(1), U, V, A );
        const Real frobA = El::FrobeniusNorm( A );
        if( print )
            El::Print( A, "A" );

        El::QRCtrl<Real> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = maxSteps;
        if( tol != Real(-1) )
        {
            ctrl.adaptive = true;
            ctrl.tol = tol;
        }
        El::DistPermutation PR(grid), PC(grid);
        El::DistMatrix<Scalar> Z(grid);
        El::Timer timer;
        if( El::mpi::Rank(comm) == 0 )
            timer.Start();
        El::Skeleton( A, PR, PC, Z, ctrl );
        if( El::mpi::Rank(comm) == 0 )
            timer.Stop();
        const El::Int rank = Z.Height();
        if( print )
        {
            El::DistMatrix<El::Int> PFull(grid);
            PR.ExplicitMatrix( PFull );
            El::Print( PFull, "PR" );
            PC.ExplicitMatrix( PFull );
            El::Print( PFull, "PC" );
            El::Print( Z, "Z" );
        }

        // Form the matrices of A's (hopefully) dominant rows and columns
        El::DistMatrix<Scalar> AR( A );
        PR.PermuteRows( AR );
        AR.Resize( rank, A.Width() );
        El::DistMatrix<Scalar> AC( A );
        PC.PermuteCols( AC );
        AC.Resize( A.Height(), rank );
        if( print )
        {
            El::Print( AC, "A_C" );
            El::Print( AR, "A_R" );
        }

        // Check || A - AC Z AR ||_F / || A ||_F
        El::DistMatrix<Scalar> B(grid);
        El::Gemm( El::NORMAL, El::NORMAL, Scalar(1), Z, AR, B );
        El::Gemm( El::NORMAL, El::NORMAL, Scalar(-1), AC, B, Scalar(1), A );
        const Real frobError = El::FrobeniusNorm( A );
        if( print )
            El::Print( A, "A - A_C Z A_R" );

        if( El::mpi::Rank(comm) == 0 )
        {
            El::Output("Skeleton time: ",timer.Total()," secs");
            El::Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A - A_C Z A_R ||_F / || A ||_F = ",frobError/frobA);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
