/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
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
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int numRhs = Input("--numRhs","number of random r.h.s.",100);
        const double realMean = Input("--realMean","real mean",0.); 
        const double imagMean = Input("--imagMean","imag mean",0.);
        const double stddev = Input("--stddev","standard dev.",1.);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        const Int pivotInt = Input("--pivot","pivot type",0);
        const double gamma = Input("--gamma","pivot constant",0.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        const auto pivotType = static_cast<LDLPivotType>(pivotInt);
        LDLPivotCtrl<Real> ctrl(pivotType);
        if( gamma != Real(0) )
            ctrl.gamma = gamma;

        C mean( realMean, imagMean );
        Matrix<C> A;
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
        Matrix<Int> perm;
        Matrix<C> dSub, factA( A );
        MakeTrapezoidal( LOWER, factA );
        LDL( factA, dSub, perm, conjugate, ctrl );
        if( print && mpi::WorldRank()==0 )
        {
            Print( A,     "A"     );
            Print( factA, "factA" );
            Print( dSub,  "dSub"  );
            Print( perm,  "perm" );
        }

        // Generate a random set of vectors
        Matrix<C> X;
        Uniform( X, n, numRhs );
        Matrix<C> B;
        Zeros( B, n, numRhs );
        Symm( LEFT, LOWER, C(1), A, X, C(0), B, conjugate );
        if( print && mpi::WorldRank()==0 )
        {
            Print( X, "X" );
            Print( B, "B" );
        }
        ldl::SolveAfter( factA, dSub, perm, B, conjugate );
        const Real AFrob = HermitianFrobeniusNorm( LOWER, A );
        const Real XFrob = FrobeniusNorm( X );
        Axpy( C(-1), B, X );
        const Real errFrob = FrobeniusNorm( X );
        if( print && mpi::WorldRank() == 0 )
        {
            Print( B, "XComp" );
            Print( X, "E" );
        }
        if( mpi::WorldRank() == 0 )
        {
            cout << "|| A ||_F = " << AFrob << "\n"
                 << "|| X ||_F = " << XFrob << "\n"
                 << "|| X - inv(A) A X ||_F = " << errFrob << endl;
        }

        if( conjugate )
        {
            // Compute the inertia of A now that we are done with it.
            auto inertia = Inertia( LOWER, A, ctrl );
            if( mpi::WorldRank() == 0 )
            {
                cout << "numPositive=" << inertia.numPositive << "\n"
                     << "numNegative=" << inertia.numNegative << "\n"
                     << "numZero    =" << inertia.numZero << "\n" << endl;
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
