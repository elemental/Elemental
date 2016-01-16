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

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int n = Input("--size","size of matrix to factor",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int numRhs = Input("--numRhs","number of random r.h.s.",100);
        const Real realMean = Input("--realMean","real mean",Real(0)); 
        const Real imagMean = Input("--imagMean","imag mean",Real(0));
        const Real stddev = Input("--stddev","standard dev.",Real(1));
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        const Int pivotInt = Input("--pivot","pivot type",0);
        const Real gamma = Input("--gamma","pivot constant",Real(0));
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        const auto pivotType = static_cast<LDLPivotType>(pivotInt);
        LDLPivotCtrl<Real> ctrl(pivotType);
        if( gamma != Real(0) )
            ctrl.gamma = gamma; 

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
        DistPermutation p;
        DistMatrix<C,MD,STAR> dSub;
        DistMatrix<C> factA( A );
        MakeTrapezoidal( LOWER, factA );
        LDL( factA, dSub, p, conjugate, ctrl );
        if( print )
        {
            Print( A,     "A"     );
            Print( factA, "factA" );
            Print( dSub,  "dSub"  );
 
            DistMatrix<Int> P;
            p.ExplicitMatrix( P );
            Print( P, "P" );
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
        X -= B;
        const Real errFrob = FrobeniusNorm( X );
        if( print )
        {
            Print( B, "XComp" );
            Print( X, "E" );
        }
        if( mpi::Rank() == 0 )
            Output
            ("|| A ||_F = ",AFrob,"\n",
             "|| X ||_F = ",XFrob,"\n",
             "|| X - inv(A) A X ||_F = ",errFrob);

        if( conjugate )
        {
            // Compute the inertia of A now that we are done with it.
            auto inertia = Inertia( LOWER, A, ctrl );
            if( mpi::Rank() == 0 )
                Output
                ("numPositive=",inertia.numPositive,"\n",
                 "numNegative=",inertia.numNegative,"\n",
                 "numZero    =",inertia.numZero,"\n");
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
