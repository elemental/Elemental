/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

const Int m=300, n=300;  // run SVD on m x n matrix
typedef double Real;     // real datatype is `Real'
typedef Complex<Real> C; // complex datatype `C'

int main()
{
    Environment env;
    DistMatrix<C> A;
    Uniform( A, m, n );

    Timer timer;
    DistMatrix<C> U, V;
    DistMatrix<Real,VR,STAR> s;
    if( mpi::Rank() == 0 )
        timer.Start();
    SVD( A, U, s, V );
    if( mpi::Rank() == 0 )
        timer.Stop();
    const Real twoNormA = MaxNorm( s );

    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
    const Real frobNormE = FrobeniusNorm( A );
    const Real eps = limits::Epsilon<Real>();
    const Real scaledResid = frobNormE / (Max(m,n)*eps*twoNormA);
    if( mpi::Rank() == 0 )
    {
        Output("SimpleSVD time: ",timer.Total()," secs");
        Output("||A||_2 = ",twoNormA);
        Output("||A - U Sigma V^H||_F / (max(m,n) eps ||A||_2) = ",scaledResid);
    }
}
