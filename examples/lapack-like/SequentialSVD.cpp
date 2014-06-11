/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_MATRICES_INC
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    enum TestType { FOURIER=0, HILBERT=1, IDENTITY=2, ONES=3, ONE_TWO_ONE=4,
                    UNIFORM=5, WILKINSON=6, ZEROS=7 }; 

    try 
    {
        const Int k = Input("--size","problem size",100);
        ProcessInput();
        PrintInputReport();

        Matrix<C> A, U, V;
        Matrix<Real> s;

        for( Int test=0; test<16; ++test )
        {
            Int n;
            const TestType testType = TestType(test/2);
            const bool useQR = test % 2;
            const std::string qrString = ( useQR ? "with QR:" : "with D&C:" );
            switch( testType )
            {
            case FOURIER:     
                if( mpi::WorldRank() == 0 ) 
                    std::cout << "Testing Fourier " 
                              << qrString << std::endl;
                n = k;
                Fourier( A, n ); 
                break;
            case HILBERT:     
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Hilbert " << qrString << std::endl;
                n = k;
                Hilbert( A, n ); 
                break;
            case IDENTITY:    
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Identity " << qrString << std::endl;
                n = k;
                Identity( A, n, n ); 
                break;
            case ONES:        
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Ones " << qrString << std::endl;
                n = k;
                Ones( A, n, n ); 
                break;
            case ONE_TWO_ONE: 
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing OneTwoOne " << qrString << std::endl;
                n = k;
                OneTwoOne( A, n ); 
                break;
            case UNIFORM:     
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Uniform " << qrString << std::endl;
                n = k;
                Uniform( A, n, n ); 
                break;
            case WILKINSON:   
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Wilkinson " << qrString << std::endl;
                Wilkinson( A, k ); 
                n = 2*k+1;
                break;
            case ZEROS:       
                if( mpi::WorldRank() == 0 )
                    std::cout << "Testing Zeros " << qrString << std::endl;
                n = k;
                Zeros( A, n, n ); 
                break;
            };

            // Make a copy of A and then perform the SVD
            U = A;
            SVD( U, s, V, useQR );

            const Real twoNormOfA = MaxNorm( s );
            const Real maxNormOfA = MaxNorm( A );
            const Real oneNormOfA = OneNorm( A );
            const Real infNormOfA = InfinityNorm( A );
            const Real frobNormOfA = FrobeniusNorm( A );
            const Real twoEstOfA = TwoNormEstimate( A );

            DiagonalScale( RIGHT, NORMAL, s, U );
            Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
            const Real maxNormOfE = MaxNorm( A );
            const Real oneNormOfE = OneNorm( A );
            const Real infNormOfE = InfinityNorm( A );
            const Real frobNormOfE = FrobeniusNorm( A );
            const Real epsilon = lapack::MachineEpsilon<Real>();
            const Real scaledResidual = frobNormOfE/(n*epsilon*twoNormOfA);

            if( mpi::WorldRank() == 0 )
            {
                cout << "||A||_max   = " << maxNormOfA << "\n"
                     << "||A||_1     = " << oneNormOfA << "\n"
                     << "||A||_oo    = " << infNormOfA << "\n"
                     << "||A||_F     = " << frobNormOfA << "\n"
                     << "\n"
                     << "||A||_2     = " << twoNormOfA << "\n"
                     << "||A||_2 est = " << twoEstOfA << "\n"
                     << "\n"
                     << "||A - U Sigma V^H||_max = " << maxNormOfE << "\n"
                     << "||A - U Sigma V^H||_1   = " << oneNormOfE << "\n"
                     << "||A - U Sigma V^H||_oo  = " << infNormOfE << "\n"
                     << "||A - U Sigma V^H||_F   = " << frobNormOfE << "\n"
                     << "||A - U Sigma V_H||_F / (n eps ||A||_2) = " 
                     << scaledResidual << "\n" << endl;
            }
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
