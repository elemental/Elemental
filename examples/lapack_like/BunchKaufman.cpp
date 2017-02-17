/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--size","size of matrix to factor",100);
        const El::Int nb = El::Input("--nb","algorithmic blocksize",96);
        const El::Int numRhs =
          El::Input("--numRhs","number of random r.h.s.",100);
        const double realMean = El::Input("--realMean","real mean",0.);
        const double imagMean = El::Input("--imagMean","imag mean",0.);
        const double stddev = El::Input("--stddev","standard dev.",1.);
        const bool conjugate = El::Input("--conjugate","LDL^H?",false);
        const El::Int pivotInt = El::Input("--pivot","pivot type",0);
        const double gamma = El::Input("--gamma","pivot constant",0.);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( nb );
        const auto pivotType = static_cast<El::LDLPivotType>(pivotInt);
        El::LDLPivotCtrl<double> ctrl(pivotType);
        if( gamma != 0. )
            ctrl.gamma = gamma;

        El::Complex<double> mean( realMean, imagMean );
        El::DistMatrix<El::Complex<double>> A;
        if( conjugate )
        {
            El::Wigner( A, n, mean, stddev );
        }
        else
        {
            El::Gaussian( A, n, n, mean, stddev );
            El::MakeSymmetric( El::LOWER, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        El::DistPermutation p;
        El::DistMatrix<El::Complex<double>> dSub, factA( A );
        El::MakeTrapezoidal( El::LOWER, factA );
        El::LDL( factA, dSub, p, conjugate, ctrl );
        if( print )
        {
            El::Print( A,     "A"     );
            El::Print( factA, "factA" );
            El::Print( dSub,  "dSub"  );

            El::DistMatrix<El::Int> P;
            p.ExplicitMatrix( P );
            El::Print( P, "P" );
        }

        // Generate a random set of vectors
        El::DistMatrix<El::Complex<double>> X;
        El::Uniform( X, n, numRhs );
        El::DistMatrix<El::Complex<double>> B;
        El::Zeros( B, n, numRhs );
        El::Symm
        ( El::LEFT, El::LOWER,
          El::Complex<double>(1), A, X, El::Complex<double>(0), B, conjugate );
        if( print )
        {
            El::Print( X, "X" );
            El::Print( B, "B" );
        }
        El::ldl::SolveAfter( factA, dSub, p, B, conjugate );
        const double AFrob = El::HermitianFrobeniusNorm( El::LOWER, A );
        const double XFrob = El::FrobeniusNorm( X );
        X -= B;
        const double errFrob = El::FrobeniusNorm( X );
        if( print )
        {
            El::Print( B, "XComp" );
            El::Print( X, "E" );
        }
        if( El::mpi::Rank() == 0 )
            El::Output
            ("|| A ||_F = ",AFrob,"\n",
             "|| X ||_F = ",XFrob,"\n",
             "|| X - inv(A) A X ||_F = ",errFrob);

        if( conjugate )
        {
            // Compute the inertia of A now that we are done with it.
            auto inertia = El::Inertia( El::LOWER, A, ctrl );
            if( El::mpi::Rank() == 0 )
                El::Output
                ("numPositive=",inertia.numPositive,"\n",
                 "numNegative=",inertia.numNegative,"\n",
                 "numZero    =",inertia.numZero,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
