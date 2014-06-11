/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_WIGNER_INC
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix to factor",100);
        const double realMean = Input("--realMean","real mean",0.); 
        const double imagMean = Input("--imagMean","imag mean",0.);
        const double stddev = Input("--stddev","standard dev.",1.);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        C mean( realMean, imagMean );
        DistMatrix<C> A;
        if( conjugate )
        {
            Wigner( A, n, mean, stddev );
        }
        else
        {
            Gaussian( A, n, n, mean, stddev );
            MakeSymmetric( LOWER, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        // WARNING: There is no pivoting here!
        DistMatrix<C> factA( A );
        LDL( factA, conjugate );
        auto d = factA.GetDiagonal();

        DistMatrix<C> L( factA );
        MakeTriangular( LOWER, L );
        SetDiagonal( L, C(1) );

        DistMatrix<C> LD( L );
        DiagonalScale( RIGHT, NORMAL, d, LD );
        Gemm( NORMAL, orientation, C(-1), LD, L, C(1), A );
        const Real frobNormOfError = FrobeniusNorm( A );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| A - L D L^[T/H] ||_F = " << frobNormOfError << "\n"
                      << std::endl;
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
