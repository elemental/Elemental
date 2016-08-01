/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    enum TestType { FOURIER=0, HILBERT=1, IDENTITY=2, ONES=3, ONE_TWO_ONE=4,
                    UNIFORM=5, WILKINSON=6, ZEROS=7 }; 

    try 
    {
        const Int k = Input("--size","problem size",100);
        const bool useLAPACK = Input("--useLAPACK","use LAPACK?",false);
        const bool useLAPACKQR = Input("--useLAPACKQR","use LAPACK QR?",false);
        const Int numTests = Input("--numTests","number of tests",5);
        ProcessInput();
        PrintInputReport();

        Matrix<C> A, U, V;
        Matrix<Real> s;

        SVDCtrl<Real> ctrl;
        for( Int test=0; test<numTests; ++test )
        {
            Int n=1;
            const TestType testType = TestType(test/2);
            ctrl.useLAPACK = useLAPACK;
            ctrl.bidiagSVDCtrl.qrCtrl.useLAPACK = useLAPACKQR;
            switch( testType )
            {
            case FOURIER:     
                if( mpi::Rank() == 0 ) 
                    Output("Testing Fourier");
                n = k;
                Fourier( A, n ); 
                break;
            case HILBERT:     
                if( mpi::Rank() == 0 )
                    Output("Testing Hilbert");
                n = k;
                Hilbert( A, n ); 
                break;
            case IDENTITY:    
                if( mpi::Rank() == 0 )
                    Output("Testing Identity");
                n = k;
                Identity( A, n, n ); 
                break;
            case ONES:        
                if( mpi::Rank() == 0 )
                    Output("Testing Ones");
                n = k;
                Ones( A, n, n ); 
                break;
            case ONE_TWO_ONE: 
                if( mpi::Rank() == 0 )
                    Output("Testing OneTwoOne");
                n = k;
                OneTwoOne( A, n ); 
                break;
            case UNIFORM:     
                if( mpi::Rank() == 0 )
                    Output("Testing Uniform");
                n = k;
                Uniform( A, n, n ); 
                break;
            case WILKINSON:   
                if( mpi::Rank() == 0 )
                    Output("Testing Wilkinson");
                Wilkinson( A, k ); 
                n = 2*k+1;
                break;
            case ZEROS:       
                if( mpi::Rank() == 0 )
                    Output("Testing Zeros");
                n = k;
                Zeros( A, n, n ); 
                break;
            };

            SVD( A, U, s, V, ctrl );

            const Real twoNormA = MaxNorm( s );
            const Real maxNormA = MaxNorm( A );
            const Real frobNormA = FrobeniusNorm( A );
            const Real twoEstA = TwoNormEstimate( A );

            DiagonalScale( RIGHT, NORMAL, s, U );
            Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
            const Real maxNormE = MaxNorm( A );
            const Real frobNormE = FrobeniusNorm( A );
            const Real epsilon = limits::Epsilon<Real>();
            const Real scaledResidual = frobNormE/(n*epsilon*twoNormA);

            if( mpi::Rank() == 0 )
                Output
                ("||A||_max   = ",maxNormA,"\n",
                 "||A||_F     = ",frobNormA,"\n",
                 "||A||_2     = ",twoNormA,"\n",
                 "||A||_2 est = ",twoEstA,"\n",
                 "||A - U Sigma V^H||_max = ",maxNormE,"\n",
                 "||A - U Sigma V^H||_F   = ",frobNormE,"\n",
                 "||A - U Sigma V_H||_F / (n eps ||A||_2) = ",scaledResidual);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
