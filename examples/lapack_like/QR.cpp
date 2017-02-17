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

        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::Int nb = El::Input("--nb","blocksize",96);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( nb );

        const El::Grid grid( comm );
        El::DistMatrix<Scalar> A(grid);
        El::Uniform( A, m, n );
        if( print )
            El::Print( A, "A" );
        const Real frobA = El::FrobeniusNorm( A );

        // Compute the QR decomposition of A, but do not overwrite A
        El::DistMatrix<Scalar> Q( A ), QFull( A ), R;
        El::qr::Explicit( Q, R );
        El::qr::Explicit( QFull, R, false );
        if( print )
        {
            El::Print( Q, "Q" );
            El::Print( QFull, "QFull" );
            El::Print( R, "R" );
        }

        // Check the error in the QR factorization, || A - Q R ||_F / || A ||_F
        El::DistMatrix<Scalar> E( A );
        El::Gemm( El::NORMAL, El::NORMAL, Scalar(-1), Q, R, Scalar(1), E );
        const Real frobQR = El::FrobeniusNorm( E );

        // Check the numerical orthogonality of Q, || I - Q^H Q ||_F / || A ||_F
        const El::Int k = El::Min(m,n);
        El::Identity( E, k, k );
        El::Herk( El::LOWER, El::ADJOINT, Real(-1), Q, Real(1), E );
        const Real frobOrthog = El::HermitianFrobeniusNorm( El::LOWER, E );

        if( El::mpi::Rank(comm) == 0 )
            El::Output
            ("|| A ||_F = ",frobA,"\n",
             "|| A - Q R ||_F / || A ||_F   = ",frobQR/frobA,"\n",
             "|| I - Q^H Q ||_F / || A ||_F = ",frobOrthog/frobA,"\n");
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
