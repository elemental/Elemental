/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try 
    {
        const Int m = Input("--height","height of matrix",20);
        const Int n = Input("--width","width of matrix",100);
        const Int r = Input("--rank","rank of matrix",5);
        Int gridHeight = Input("--gridHeight","grid height",0);
        const Int pnp = Input("--parnp","number of particitaed processes",2);
        const Int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const double tol = Input("--tol","tolerance for ID",-1.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( commSize <= pnp )
        {
            if( commRank == 0 )
                std::cout << "Do not have enough processes to be viewers." 
                          << std::endl;
            Finalize();
            return 0;
        }

        if( gridHeight == 0 )
            gridHeight = Grid::FindFactor( pnp );

        mpi::Group group;
        mpi::CommGroup( comm, group );
        std::vector<Int> ranks(pnp);
        mpi::Group pgroup;
        for( int q=0; q<pnp; ++q )
            ranks[q] = q;
        mpi::Incl( group, pnp, ranks.data(), pgroup );
        Grid g( comm, pgroup, gridHeight );

        Log( "Gemm, A = UV^*" );
        DistMatrix<C> U(g), V(g), A(g);
        Uniform( U, m, r );
        Uniform( V, n, r );
        A.Resize( m, n );
        if( commRank < pnp )
            Gemm( NORMAL, ADJOINT, C(1), U, V, A );
        
        // Invalid root error message
        Log( "Test of auto AA(A) with viewers" );
        auto AA(A);

        Log( "Forbenious Norm of A" );
        Real frobA = Real(0.0);
        if( commRank < pnp )
            frobA = FrobeniusNorm( A );
        if( print && commRank < pnp )
            Print( A, "A", LogOS() );

        Log( "ID of A" );
        DistMatrix<Int,VR,STAR> p(g);
        DistMatrix<C> Z(g);
        QRCtrl<double> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = maxSteps;
        if( tol != -1. )
        {
            ctrl.adaptive = true;
            ctrl.tol = tol;
        }
        if( commRank < pnp )
            ID( A, p, Z, ctrl );
        const Int rank = Z.Height();
        if( print )
        {
            Print( p, "p" );
            Print( Z, "Z" );
        }

        /* We do not currently care about the accuracy of ID
        // Pivot A and form the matrix of its (hopefully) dominant columns
        Log( "PermuteCol of A" );
        if( commRank < pnp )
            InversePermuteCols( A, p );
        Log( "hatA = A" );
        auto hatA( A );
        hatA.Resize( m, rank );
        if( print && commRank < pnp )
        {
            Print( A, "A P", LogOS() );
            Print( hatA, "\\hat{A}", LogOS() );
        }

        // Check || A P - \hat{A} [I, Z] ||_F / || A ||_F
        Log( "Prtition A = [AL AR]" );
        DistMatrix<C> AL(g), AR(g);
        if( commRank < pnp )
            PartitionRight( A, AL, AR, rank );
        if( commRank < pnp )
            Zero( AL );
        if( commRank < pnp )
        {
            DistMatrix<C,MC,STAR> hatA_MC_STAR(g);
            DistMatrix<C,STAR,MR> Z_STAR_MR(g);
            hatA_MC_STAR.AlignWith( AR );
            Z_STAR_MR.AlignWith( AR );
            hatA_MC_STAR = hatA;
            Z_STAR_MR = Z;
            LocalGemm
            ( NORMAL, NORMAL, C(-1), hatA_MC_STAR, Z_STAR_MR, C(1), AR );
        }
        Real frobError = Real(0.0);
        if( commRank < pnp )
            frobError = FrobeniusNorm( A );
        if( print && commRank < pnp )
            Print( A, "A P - \\hat{A} [I, Z]", LogOS() );

        Log( "Show error" );
        if( commRank < pnp )
            LogOS() << "|| A ||_F = " << frobA << "\n\n"               
                    << "|| A P - \\hat{A} [I, Z] ||_F / || A ||_F = " 
                    << frobError/frobA << "\n" << std::endl;
                    */
        Log( "Finished" );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
