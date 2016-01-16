/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> F;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",20);
        const Int n = Input("--width","width of matrix",100);
        const Int r = Input("--rank","rank of matrix",5);
        const Int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const Real tol = Input("--tol","tolerance for ID",Real(-1));
        const bool print = Input("--print","print matrices?",false);
        const bool smallestFirst =
          Input("--smallestFirst","smallest norm first?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        DistMatrix<F> U(g), V(g), A(g);
        Uniform( U, m, r );
        Uniform( V, n, r );
        Gemm( NORMAL, ADJOINT, F(1), U, V, A );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            Print( A, "A" );

        DistPermutation Omega(g);
        DistMatrix<F,STAR,VR> Z(g);
        QRCtrl<Real> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = maxSteps;
        if( tol != -1. )
        {
            ctrl.adaptive = true;
            ctrl.tol = tol;
        }
        ctrl.smallestFirst = smallestFirst;
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        ID( A, Omega, Z, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        const Int rank = Z.Height();
        if( print )
        {
            DistMatrix<Int> OmegaFull(g);
            Omega.ExplicitMatrix( OmegaFull );
            Print( OmegaFull, "Omega" );
            Print( Z, "Z" );
        }

        // Pivot A and form the matrix of its (hopefully) dominant columns
        Timer permTimer( "permTimer" );
        permTimer.Start();
        Omega.PermuteCols( A );
        permTimer.Stop();

        auto hatA( A );
        hatA.Resize( m, rank );
        if( print )
        {
            Print( A, "A Omega^T" );
            Print( hatA, "\\hat{A}" );
        }

        // Check || A Omega^T - \hat{A} [I, Z] ||_F / || A ||_F
        DistMatrix<F> AL(g), AR(g);
        PartitionRight( A, AL, AR, rank );
        Zero( AL );
        {
            DistMatrix<F,MC,STAR> hatA_MC_STAR(g);
            DistMatrix<F,STAR,MR> Z_STAR_MR(g);
            hatA_MC_STAR.AlignWith( AR );
            Z_STAR_MR.AlignWith( AR );
            hatA_MC_STAR = hatA;
            Z_STAR_MR = Z;
            LocalGemm
            ( NORMAL, NORMAL, F(-1), hatA_MC_STAR, Z_STAR_MR, F(1), AR );
        }
        const Real frobError = FrobeniusNorm( A );
        if( print )
            Print( A, "A Omega^T - \\hat{A} [I, Z]" );

        if( mpi::Rank() == 0 )
        {
            Output("  ID time: ",timer.Total()," secs");
            Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A Omega^T - \\hat{A} [I, Z] ||_F / || A ||_F = ",
             frobError/frobA);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
