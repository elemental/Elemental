/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","height of matrix",20);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::Int r = El::Input("--rank","rank of matrix",5);
        const El::Int maxSteps =
          El::Input("--maxSteps","max # of steps of QR",10);
        const Real tol = El::Input("--tol","tolerance for ID",Real(-1));
        const bool print = El::Input("--print","print matrices?",false);
        const bool smallestFirst =
          El::Input("--smallestFirst","smallest norm first?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::mpi::Comm comm = El::mpi::COMM_WORLD;

        const El::Grid& grid = El::Grid::Default();
        El::DistMatrix<Scalar> U(grid), V(grid), A(grid);
        El::Uniform( U, m, r );
        El::Uniform( V, n, r );
        El::Gemm( El::NORMAL, El::ADJOINT, Scalar(1), U, V, A );
        const Real frobA = El::FrobeniusNorm( A );
        if( print )
            El::Print( A, "A" );

        El::DistPermutation Omega(grid);
        El::DistMatrix<Scalar,El::STAR,El::VR> Z(grid);
        El::QRCtrl<Real> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = maxSteps;
        if( tol != -1. )
        {
            ctrl.adaptive = true;
            ctrl.tol = tol;
        }
        ctrl.smallestFirst = smallestFirst;
        El::Timer timer;
        if( El::mpi::Rank(comm) == 0 )
            timer.Start();
        El::ID( A, Omega, Z, ctrl );
        if( El::mpi::Rank(comm) == 0 )
            timer.Stop();
        const El::Int rank = Z.Height();
        if( print )
        {
            El::DistMatrix<El::Int> OmegaFull(grid);
            Omega.ExplicitMatrix( OmegaFull );
            El::Print( OmegaFull, "Omega" );
            El::Print( Z, "Z" );
        }

        // Pivot A and form the matrix of its (hopefully) dominant columns
        El::Timer permTimer("permTimer");
        permTimer.Start();
        Omega.PermuteCols( A );
        permTimer.Stop();

        auto hatA( A );
        hatA.Resize( m, rank );
        if( print )
        {
            El::Print( A, "A Omega^T" );
            El::Print( hatA, "\\hat{A}" );
        }

        // Check || A Omega^T - \hat{A} [I, Z] ||_F / || A ||_F
        El::DistMatrix<Scalar> AL(grid), AR(grid);
        El::PartitionRight( A, AL, AR, rank );
        El::Zero( AL );
        {
            El::DistMatrix<Scalar,El::MC,El::STAR> hatA_MC_STAR(grid);
            El::DistMatrix<Scalar,El::STAR,El::MR> Z_STAR_MR(grid);
            hatA_MC_STAR.AlignWith( AR );
            Z_STAR_MR.AlignWith( AR );
            hatA_MC_STAR = hatA;
            Z_STAR_MR = Z;
            El::LocalGemm
            ( El::NORMAL, El::NORMAL,
              Scalar(-1), hatA_MC_STAR, Z_STAR_MR, Scalar(1), AR );
        }
        const Real frobError = El::FrobeniusNorm( A );
        if( print )
            El::Print( A, "A Omega^T - \\hat{A} [I, Z]" );

        if( El::mpi::Rank(comm) == 0 )
        {
            El::Output("  ID time: ",timer.Total()," secs");
            El::Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A Omega^T - \\hat{A} [I, Z] ||_F / || A ||_F = ",
             frobError/frobA);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
