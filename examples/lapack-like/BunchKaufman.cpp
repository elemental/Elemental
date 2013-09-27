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
#include "elemental/blas-like/level1/MakeSymmetric.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level1/Transpose.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors/Util.hpp"
#include "elemental/lapack-like/LDL.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/HermitianUniformSpectrum.hpp"
using namespace std;
using namespace elem;

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
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double realMean = Input("--realMean","real mean",0.); 
        const double imagMean = Input("--imagMean","imag mean",0.);
        const double stddev = Input("--stddev","standard dev.",1.);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        C mean( realMean, imagMean );
        DistMatrix<C> A;
        if( conjugate )
        {
            // Wigner yields an excessively high (and inaccurate) determinant
            //Wigner( A, n, mean, stddev );
            HermitianUniformSpectrum( A, n, 1., 2. );
        }
        else
        {
            Gaussian( A, n, n, mean, stddev );
            MakeSymmetric( LOWER, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        DistMatrix<Int,VC,STAR> p;
        DistMatrix<C> factA( A );
        MakeTriangular( LOWER, factA );
        if( conjugate )
            LDLH( factA, p );
        else
            LDLT( factA, p );
        Print( A,     "A"     );
        Print( factA, "factA" );
        Print( p,     "p"     );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
