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
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        ProcessInput();
        PrintInputReport();

        Matrix<C> A;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );

        // Compute the QR decomposition of A, but do not overwrite A
        Matrix<C> Q( A ), R;
        qr::Explicit( Q, R );

        // Check the error in the QR factorization, || A - Q R ||_F / || A ||_F
        Matrix<C> E( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, R, C(1), E );
        const Real frobQR = FrobeniusNorm( E );

        // Check the numerical orthogonality of Q, || I - Q^H Q ||_F / || A ||_F
        const Int k = Min(m,n);
        Identity( E, k, k );
        Herk( LOWER, ADJOINT, Real(-1), Q, Real(1), E );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, E );

        if( mpi::Rank() == 0 )
            Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A - Q R ||_F / || A ||_F   = ",frobQR/frobA,"\n",
             "|| I - Q^H Q ||_F / || A ||_F = ",frobOrthog/frobA,"\n");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
