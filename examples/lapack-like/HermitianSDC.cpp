/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HermitianEig.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Wigner.hpp"
#include "elemental/io/Display.hpp"
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int n = Input("--size","height of matrix",100);
        const Int cutoff = Input("--cutoff","cutoff for QR alg.",256);
        const Int maxInnerIts = Input("--maxInnerIts","maximum RURV its",1);
        const Int maxOuterIts = Input("--maxOuterIts","maximum it's/split",10);
        const Real relTol = Input("--relTol","rel. tol.",Real(0));
        const bool display = Input("--display","display matrices?",false);
        ProcessInput();
        PrintInputReport();

        const Grid& g = DefaultGrid();
        auto A = Wigner<C>( g, n );
        const Real frobA = FrobeniusNorm( A );

        // Attempt to compute the spectral decomposition of A, 
        // but do not overwrite A
        DistMatrix<C> D( A ), Q(g);
        hermitian_eig::SDC
        ( LOWER, D, Q, cutoff, maxInnerIts, maxOuterIts, relTol );
        MakeTriangular( LOWER, D ); 

        if( display )
        {
            Display( A, "A" );
            Display( Q, "Q" );
            Display( D.GetRealPartOfDiagonal(), "w" );
        }

        DistMatrix<C> G(g);
        Gemm( NORMAL, NORMAL, C(1), Q, D, G );
        Gemm( NORMAL, ADJOINT, C(-1), G, Q, C(1), A );
        MakeTrapezoidal( LOWER, D, -1 );
        const Real frobOffD = HermitianFrobeniusNorm( LOWER, D );
        const Real frobE = FrobeniusNorm( A ); 
        if( mpi::WorldRank() == 0 )
        {
            std::cout << " || A - Q D Q^H ||_F / || A ||_F = " << frobE/frobA 
                      << "\n"
                      << " || off(D) ||_F      / || A ||_F = " << frobOffD/frobA
                      << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
