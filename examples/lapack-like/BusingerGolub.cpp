/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/UpdateDiagonal.hpp"
#include "elemental/lapack-like/ApplyColumnPivots.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
#include "elemental/io.hpp"
using namespace std;
using namespace elem;

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
        const bool alwaysRecompute = Input("--always","no norm updates?",false);
        const bool blockedUnpiv = 
            Input("--blockUnpiv","blocked unpivoted QR?",false);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        auto A = Uniform<C>( g, m, n );
        const Real frobA = FrobeniusNorm( A );
        if( display )
            Display( A, "A" );
        if( print )
            Print( A, "A" );

        // Compute the pivoted QR decomposition of A, but do not overwrite A
        DistMatrix<C> QRPiv( A );
        DistMatrix<C,MD,STAR> tPiv(g);
        DistMatrix<Int,VR,STAR> p(g);
        qr::BusingerGolub( QRPiv, tPiv, p, alwaysRecompute );
        if( display )
        {
            Display( QRPiv, "QRPiv" );
            Display( tPiv, "tPiv" );
            Display( p, "p" );
        }
        if( print )
        {
            Print( QRPiv, "QRPiv" );
            Print( tPiv, "tPiv" );
            Print( p, "p" );
        }

        // Compute the standard QR decomposition of A
        auto QRNoPiv( A );
        DistMatrix<C,MD,STAR> tNoPiv(g);
        if( blockedUnpiv )
            QR( QRNoPiv, tNoPiv );
        else
            qr::PanelHouseholder( QRNoPiv, tNoPiv );
        if( display )
        {
            Display( QRNoPiv, "QRNoPiv" );
            Display( tNoPiv, "tNoPiv" );
        }
        if( print )
        {
            Print( QRNoPiv, "QRNoPiv" );
            Print( tNoPiv, "tNoPiv" );
        }

        // Check the error in the pivoted QR factorization, 
        // || A P - Q R ||_F / || A ||_F
        auto E( QRPiv );
        MakeTriangular( UPPER, E );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, E );
        ApplyInverseColumnPivots( E, p ); 
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
        qr::ApplyQ( LEFT, NORMAL, QRNoPiv, tNoPiv, E );
        Axpy( C(-1), A, E );
        const Real frobQRNoPiv = FrobeniusNorm( E );
        if( display )
            Display( E, "A - Q R" );
        if( print )
            Print( E, "A - Q R" );

        // Check orthogonality of pivoted Q, || I - Q^H Q ||_F / || A ||_F
        Identity( E, m, n );
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, E );
        qr::ApplyQ( LEFT, ADJOINT, QRPiv, tPiv, E );
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
        qr::ApplyQ( LEFT, NORMAL, QRPiv, tPiv, E );
        qr::ApplyQ( LEFT, ADJOINT, QRPiv, tPiv, E );
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
