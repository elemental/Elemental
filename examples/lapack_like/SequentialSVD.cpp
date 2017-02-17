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

    enum TestType { FOURIER=0, HILBERT=1, IDENTITY=2, ONES=3, ONE_TWO_ONE=4,
                    UNIFORM=5, WILKINSON=6, ZEROS=7 };

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int k = El::Input("--size","problem size",100);
        const bool useLAPACK = El::Input("--useLAPACK","use LAPACK?",false);
        const bool useLAPACKQR =
          El::Input("--useLAPACKQR","use LAPACK QR?",false);
        const El::Int numTests = El::Input("--numTests","number of tests",5);
        El::ProcessInput();
        El::PrintInputReport();

        El::Matrix<Scalar> A, U, V;
        El::Matrix<Real> s;

        El::SVDCtrl<Real> ctrl;
        for( El::Int test=0; test<numTests; ++test )
        {
            El::Int n=1;
            const TestType testType = TestType(test/2);
            ctrl.useLAPACK = useLAPACK;
            ctrl.bidiagSVDCtrl.qrCtrl.useLAPACK = useLAPACKQR;
            switch( testType )
            {
            case FOURIER:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Fourier");
                n = k;
                El::Fourier( A, n );
                break;
            case HILBERT:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Hilbert");
                n = k;
                El::Hilbert( A, n );
                break;
            case IDENTITY:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Identity");
                n = k;
                El::Identity( A, n, n );
                break;
            case ONES:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Ones");
                n = k;
                El::Ones( A, n, n );
                break;
            case ONE_TWO_ONE:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing OneTwoOne");
                n = k;
                El::OneTwoOne( A, n );
                break;
            case UNIFORM:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Uniform");
                n = k;
                El::Uniform( A, n, n );
                break;
            case WILKINSON:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Wilkinson");
                El::Wilkinson( A, k );
                n = 2*k+1;
                break;
            case ZEROS:
                if( El::mpi::Rank(comm) == 0 )
                    El::Output("Testing Zeros");
                n = k;
                El::Zeros( A, n, n );
                break;
            };

            El::SVD( A, U, s, V, ctrl );

            const Real twoNormA = El::MaxNorm( s );
            const Real maxNormA = El::MaxNorm( A );
            const Real frobNormA = El::FrobeniusNorm( A );
            const Real twoEstA = El::TwoNormEstimate( A );

            El::DiagonalScale( El::RIGHT, El::NORMAL, s, U );
            El::Gemm( El::NORMAL, El::ADJOINT, Scalar(-1), U, V, Scalar(1), A );
            const Real maxNormE = El::MaxNorm( A );
            const Real frobNormE = El::FrobeniusNorm( A );
            const Real epsilon = El::limits::Epsilon<Real>();
            const Real scaledResidual = frobNormE/(n*epsilon*twoNormA);

            if( El::mpi::Rank(comm) == 0 )
                El::Output
                ("||A||_max   = ",maxNormA,"\n",
                 "||A||_F     = ",frobNormA,"\n",
                 "||A||_2     = ",twoNormA,"\n",
                 "||A||_2 est = ",twoEstA,"\n",
                 "||A - U Sigma V^H||_max = ",maxNormE,"\n",
                 "||A - U Sigma V^H||_F   = ",frobNormE,"\n",
                 "||A - U Sigma V_H||_F / (n eps ||A||_2) = ",scaledResidual);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
