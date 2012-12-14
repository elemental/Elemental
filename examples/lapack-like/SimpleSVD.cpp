/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace elem;

const int m=300, n=300;  // run SVD on m x n matrix
typedef double R;        // real datatype is `R'
typedef Complex<R> C;    // complex datatype `C'

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    DistMatrix<C> A;
    Uniform( m, n, A );

    DistMatrix<C> U, V;
    DistMatrix<R,VR,STAR> s;
    U = A;
    SVD( U, s, V );
    const R twoNormOfA = Norm( s, MAX_NORM );

    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
    const R frobNormOfE = Norm( A, FROBENIUS_NORM );
    const R eps = lapack::MachineEpsilon<R>();
    const R scaledResidual = frobNormOfE / (std::max(m,n)*eps*twoNormOfA);
    if( mpi::WorldRank() == 0 )
    {
        std::cout << "||A||_2 = " << twoNormOfA << "\n"
                  << "||A - U Sigma V^H||_F / (max(m,n) eps ||A||_2) = " 
                  << scaledResidual << "\n" << std::endl;
    }

    Finalize();
    return 0;
}
