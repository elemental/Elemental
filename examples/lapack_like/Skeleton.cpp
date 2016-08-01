/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",20);
        const Int n = Input("--width","width of matrix",100);
        Int targetRank = Input("--rank","rank of matrix",5);
        Int maxSteps = Input("--maxSteps","max # of steps of QR",10);
        const Real tol = Input("--tol","tolerance for ID",Real(-1));
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Int minDim = Min( m, n );
        targetRank = Min( targetRank, minDim );
        maxSteps = Min( maxSteps, targetRank );

        DistMatrix<C> U, V;
        Uniform( U, m, targetRank );
        Uniform( V, n, targetRank );
        DistMatrix<C> A;
        Gemm( NORMAL, ADJOINT, C(1), U, V, A );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            Print( A, "A" );

        const Grid& g = A.Grid();
        QRCtrl<Real> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = maxSteps;
        if( tol != Real(-1) )
        {
            ctrl.adaptive = true;
            ctrl.tol = tol;
        }
        DistPermutation PR(g), PC(g);
        DistMatrix<C> Z(g);
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        Skeleton( A, PR, PC, Z, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        const Int rank = Z.Height();
        if( print )
        {
            DistMatrix<Int> PFull(g);
            PR.ExplicitMatrix( PFull );
            Print( PFull, "PR" );
            PC.ExplicitMatrix( PFull );
            Print( PFull, "PC" );
            Print( Z, "Z" );
        }

        // Form the matrices of A's (hopefully) dominant rows and columns
        DistMatrix<C> AR( A );
        PR.PermuteRows( AR );
        AR.Resize( rank, A.Width() );
        DistMatrix<C> AC( A );
        PC.PermuteCols( AC );
        AC.Resize( A.Height(), rank );
        if( print )
        {
            Print( AC, "A_C" );
            Print( AR, "A_R" );
        }

        // Check || A - AC Z AR ||_F / || A ||_F
        DistMatrix<C> B(g);
        Gemm( NORMAL, NORMAL, C(1), Z, AR, B );
        Gemm( NORMAL, NORMAL, C(-1), AC, B, C(1), A );
        const Real frobError = FrobeniusNorm( A );
        if( print )
            Print( A, "A - A_C Z A_R" );

        if( mpi::Rank() == 0 )
        {
            Output("Skeleton time: ",timer.Total()," secs");
            Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A - A_C Z A_R ||_F / || A ||_F = ",frobError/frobA);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
