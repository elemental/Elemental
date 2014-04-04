/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_AXPY_INC
#include ELEM_DIAGONALSCALE_INC
#include ELEM_MAKESYMMETRIC_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_TRANSPOSE_INC
#include ELEM_SYMM_INC

#include ELEM_APPLYPACKEDREFLECTORS_INC
#include ELEM_LDL_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INERTIA_INC

#include ELEM_UNIFORM_INC
#include ELEM_WIGNER_INC
#include ELEM_ZEROS_INC
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
        const Int numRhs = Input("--numRhs","number of random r.h.s.",100);
        const double realMean = Input("--realMean","real mean",0.); 
        const double imagMean = Input("--imagMean","imag mean",0.);
        const double stddev = Input("--stddev","standard dev.",1.);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        const Int pivotInt = Input("--pivot","pivot type",0);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        const auto pivotType = static_cast<LDLPivotType>(pivotInt);

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
        DistMatrix<Int,VC,STAR> p;
        DistMatrix<C,MD,STAR> dSub;
        DistMatrix<C> factA( A );
        MakeTriangular( LOWER, factA );
        if( conjugate )
            LDLH( factA, dSub, p, pivotType );
        else
            LDLT( factA, dSub, p, pivotType );
        if( print )
        {
            Print( A,     "A"     );
            Print( factA, "factA" );
            Print( dSub,  "dSub"  );
            Print( p,     "p"     );
        }

        // Generate a random set of vectors
        DistMatrix<C> X;
        Uniform( X, n, numRhs );
        DistMatrix<C> B;
        Zeros( B, n, numRhs );
        Symm( LEFT, LOWER, C(1), A, X, C(0), B, conjugate );
        if( print )
        {
            Print( X, "X" );
            Print( B, "B" );
        }
        ldl::SolveAfter( factA, dSub, p, B, conjugate );
        const Real AFrob = HermitianFrobeniusNorm( LOWER, A );
        const Real XFrob = FrobeniusNorm( X );
        Axpy( C(-1), B, X );
        const Real errFrob = FrobeniusNorm( X );
        if( print )
        {
            Print( B, "XComp" );
            Print( X, "E" );
        }
        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| A ||_F = " << AFrob << "\n"
                      << "|| X ||_F = " << XFrob << "\n"
                      << "|| X - inv(A) A X ||_F = " << errFrob << std::endl;
        }

        if( conjugate )
        {
            // Compute the inertia of A now that we are done with it.
            auto inertia = Inertia( LOWER, A, pivotType );
            if( mpi::WorldRank() == 0 )
            {
                std::cout << "numPositive=" << inertia.numPositive << "\n"
                          << "numNegative=" << inertia.numNegative << "\n"
                          << "numZero    =" << inertia.numZero << "\n"
                          << std::endl;
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
