/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors/Util.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const int n = Input("--size","size of matrix to factor",100);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Grid g( mpi::COMM_WORLD );
        DistMatrix<C> A( g );
        if( conjugate )
        {
            HermitianUniformSpectrum( A, n, -30, -20 );
        }
        else
        {
            Uniform( A, n, n );
            DistMatrix<C> ATrans( g );
            Transpose( A, ATrans );
            Axpy( C(1), ATrans, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        // WARNING: There is no pivoting here!
        DistMatrix<C> factA( A );
        DistMatrix<C,MC,STAR> d( g );
        if( conjugate )
            LDLH( factA, d );
        else
            LDLT( factA, d );

        DistMatrix<C> L( factA );
        MakeTriangular( LOWER, L );
        SetDiagonal( LEFT, 0, L, C(1) );

        DistMatrix<C> LD( L );
        DiagonalScale( RIGHT, NORMAL, d, LD );
        Gemm( NORMAL, orientation, C(-1), LD, L, C(1), A );
        const R frobNormOfError = FrobeniusNorm( A );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| A - L D L^[T/H] ||_F = " << frobNormOfError << "\n"
                      << std::endl;
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
