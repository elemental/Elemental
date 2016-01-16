/*
   Copyright (c) 2009-2016, Jack Poulson
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
    Environment env( argc, argv );

    try 
    {
        const Int n = Input("--size","width of matrix",100);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices?",false);
        const bool smallestFirst =
          Input("--smallestFirst","smallest norm first?",false);
        ProcessInput();
        PrintInputReport();
        const Int m = n;

        QRCtrl<Real> ctrl;
        ctrl.smallestFirst = smallestFirst;

        DistMatrix<C> A;
        GKS( A, n );
        const Real frobA = FrobeniusNorm( A );
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        // Compute the pivoted QR decomposition of A, but do not overwrite A
        auto QRPiv( A );
        DistMatrix<C,MD,STAR> tPiv;
        DistMatrix<Real,MD,STAR> dPiv;
        DistPermutation Omega;
        Timer PQRtimer;
        if( mpi::Rank() == 0 )
            PQRtimer.Start();
        QR( QRPiv, tPiv, dPiv, Omega );
        if( mpi::Rank() == 0 )
            PQRtimer.Stop();
        if( display )
        {
            Display( QRPiv, "QRPiv" );
            Display( tPiv, "tPiv" );
            Display( dPiv, "dPiv" );

            DistMatrix<Int> OmegaFull;
            Omega.ExplicitMatrix( OmegaFull );
            Display( OmegaFull, "Omega" );
        }
        if( print )
        {
            Print( QRPiv, "QRPiv" );
            Print( tPiv, "tPiv" );
            Print( dPiv, "dPiv" );

            DistMatrix<Int> OmegaFull;
            Omega.ExplicitMatrix( OmegaFull );
            Print( OmegaFull, "Omega" );
        }

        // Compute the standard QR decomposition of A
        auto QRNoPiv( A );
        DistMatrix<C,MD,STAR> tNoPiv;
        DistMatrix<Real,MD,STAR> dNoPiv;
        Timer QRtimer;
        if( mpi::Rank() == 0 )
            QRtimer.Start();
        QR( QRNoPiv, tNoPiv, dNoPiv );
        if( mpi::Rank() == 0 )
            QRtimer.Start();
        if( display )
        {
            Display( QRNoPiv, "QRNoPiv" );
            Display( tNoPiv, "tNoPiv" );
            Display( dNoPiv, "dNoPiv" );
        }
        if( print )
        {
            Print( QRNoPiv, "QRNoPiv" );
            Print( tNoPiv, "tNoPiv" );
            Print( dNoPiv, "dNoPiv" );
        }

        // Check the error in the pivoted QR factorization, 
        // || A Omega^T - Q R ||_F / || A ||_F
        auto E( QRPiv );
        MakeTrapezoidal( UPPER, E );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, dPiv, E );
        Omega.InversePermuteCols( E );
        E -= A;
        const Real frobQRPiv = FrobeniusNorm( E );
        if( display )
            Display( E, "A P - Q R" );
        if( print )
            Print( E, "A P - Q R" );

        // Check the error in the standard QR factorization, 
        // || A - Q R ||_F / || A ||_F
        E = QRNoPiv;
        MakeTrapezoidal( UPPER, E );
        qr::ApplyQ( LEFT, NORMAL, QRNoPiv, tNoPiv, dNoPiv, E );
        E -= A;
        const Real frobQRNoPiv = FrobeniusNorm( E );
        if( display )
            Display( E, "A - Q R" );
        if( print )
            Print( E, "A - Q R" );

        // Check orthogonality of pivoted Q, || I - Q^H Q ||_F / || A ||_F
        Identity( E, m, n );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, dPiv, E );
        qr::ApplyQ( LEFT, ADJOINT, QRPiv, tPiv, dPiv, E );
        const Int k = Min(m,n);
        auto EUpper = View( E, 0, 0, k, k );
        ShiftDiagonal( EUpper, C(-1) );
        const Real frobOrthogPiv = FrobeniusNorm( EUpper ); 
        if( display )
            Display( E, "pivoted I - Q^H Q" );
        if( print )
            Print( E, "pivoted I - Q^H Q" );

        // Check orthogonality of unpivoted Q, || I - Q^H Q ||_F / || A ||_F
        Identity( E, m, n );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, dPiv, E );
        qr::ApplyQ( LEFT, ADJOINT, QRPiv, tPiv, dPiv, E );
        EUpper = View( E, 0, 0, k, k );
        ShiftDiagonal( EUpper, C(-1) );
        const Real frobOrthogNoPiv = FrobeniusNorm( EUpper ); 
        if( display )
            Display( E, "unpivoted I - Q^H Q" );
        if( print )
            Print( E, "unpivoted I - Q^H Q" );

        if( mpi::Rank() == 0 )
        {
            Output("Pivot QR time: ",PQRtimer.Total()," secs\n",
                   "      QR time: ",QRtimer.Total()," secs");
            Output
            ("|| A ||_F = ",frobA,"\n\n",
             "With pivoting: \n", 
             "    || A P - Q R ||_F / || A ||_F = ",frobQRPiv/frobA,"\n", 
             "    || I - Q^H Q ||_F / || A ||_F = ",frobOrthogPiv/frobA,"\n\n",
             "Without pivoting: \n",
             "    || A - Q R ||_F / || A ||_F = ",frobQRNoPiv/frobA,"\n",
             "    || I - Q^H Q ||_F / || A ||_F = ",frobOrthogNoPiv/frobA,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
