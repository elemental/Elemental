/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

#include EL_IDENTITY_INC

using namespace std;
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    const Int worldRank = mpi::WorldRank();

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        // Compute the pivoted QR decomposition of A, but do not overwrite A
        auto QRPiv( A );
        DistMatrix<C,MD,STAR> tPiv;
        DistMatrix<Real,MD,STAR> dPiv;
        DistMatrix<Int,VR,STAR> perm;
        QR( QRPiv, tPiv, dPiv, perm );
        if( display )
        {
            Display( QRPiv, "QRPiv" );
            Display( tPiv, "tPiv" );
            Display( dPiv, "dPiv" );
            Display( perm, "perm" );
        }
        if( print )
        {
            Print( QRPiv, "QRPiv" );
            Print( tPiv, "tPiv" );
            Print( dPiv, "dPiv" );
            Print( perm, "perm" );
        }

        // Compute the standard QR decomposition of A
        auto QRNoPiv( A );
        DistMatrix<C,MD,STAR> tNoPiv;
        DistMatrix<Real,MD,STAR> dNoPiv;
        QR( QRNoPiv, tNoPiv, dNoPiv );
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
        // || A P - Q R ||_F / || A ||_F
        auto E( QRPiv );
        MakeTriangular( UPPER, E );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, dPiv, E );
        InversePermuteCols( E, perm );
        Axpy( C(-1), A, E );
        const Real frobQRPiv = FrobeniusNorm( E );
        if( display )
            Display( E, "A P - Q R" );
        if( print )
            Print( E, "A P - Q R" );

        // Check the error in the standard QR factorization, 
        // || A - Q R ||_F / || A ||_F
        E = QRNoPiv;
        MakeTriangular( UPPER, E );
        qr::ApplyQ( LEFT, NORMAL, QRNoPiv, tNoPiv, dNoPiv, E );
        Axpy( C(-1), A, E );
        const Real frobQRNoPiv = FrobeniusNorm( E );
        if( display )
            Display( E, "A - Q R" );
        if( print )
            Print( E, "A - Q R" );

        // Check orthogonality of pivoted Q, || I - Q^H Q ||_F / || A ||_F
        Identity( E, m, n );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, dPiv, E );
        qr::ApplyQ( LEFT, ADJOINT, QRPiv, tPiv, dPiv, E );
        const Int k = std::min(m,n);
        auto EUpper = View( E, 0, 0, k, k );
        UpdateDiagonal( EUpper, C(-1) );
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
        UpdateDiagonal( EUpper, C(-1) );
        const Real frobOrthogNoPiv = FrobeniusNorm( EUpper ); 
        if( display )
            Display( E, "unpivoted I - Q^H Q" );
        if( print )
            Print( E, "unpivoted I - Q^H Q" );

        if( worldRank == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n\n"
                      << "With pivoting: \n" 
                      << "    || A P - Q R ||_F / || A ||_F = " 
                      << frobQRPiv/frobA << "\n"
                      << "    || I - Q^H Q ||_F / || A ||_F = "
                      << frobOrthogPiv/frobA << "\n\n"
                      << "Without pivoting: \n"
                      << "    || A - Q R ||_F / || A ||_F = "
                      << frobQRNoPiv/frobA << "\n"
                      << "    || I - Q^H Q ||_F / || A ||_F = "
                      << frobOrthogNoPiv/frobA << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
