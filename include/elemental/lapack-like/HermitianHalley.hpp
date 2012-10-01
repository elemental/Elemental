/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {

//
// A simplification of QDWH.hpp, which is loosely based on Yuji Nakatsukasa's 
// implementation of a QR-based dynamically weighted Halley iteration for the 
// polar decomposition. 
//     http://www.mathworks.com/matlabcentral/fileexchange/36830
//

template<typename F>
int HermitianHalley
( UpperOrLower uplo, Matrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("HermitianHalley");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Height must equal width");
    typedef typename Base<F>::type R;
    const int height = A.Height();
    const R oneHalf = R(1)/R(2);
    const R oneThird = R(1)/R(3);

    const R epsilon = lapack::MachineEpsilon<R>();
    const R tol = 5*epsilon;
    const R cubeRootTol = Pow(tol,oneThird);
    const R a = 3;
    const R b = 1;
    const R c = 3;

    // Form the first iterate
    Scale( 1/upperBound, A );

    int numIts=0;
    R frobNormADiff;
    Matrix<F> ALast;
    Matrix<F> Q( 2*height, height );
    Matrix<F> QT, QB;
    PartitionDown( Q, QT,
                      QB, height );
    Matrix<F> C;
    Matrix<F> ATemp;
    do
    {
        if( numIts > 100 )
            throw std::runtime_error("Halley iteration did not converge");
        ++numIts;
        ALast = A;

        // TODO: Come up with a test for when we can use the Cholesky approach
        if( true )
        {
            //
            // The standard QR-based algorithm
            //
            MakeHermitian( uplo, A );
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            ExplicitQR( Q );
            Trrk
            ( uplo, NORMAL, ADJOINT, F(a-b/c)/Sqrt(c), QT, QB, F(b/c), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            MakeHermitian( uplo, A );
            Identity( height, height, C );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( b/c, A );
            Axpy( a-b/c, ATemp, A );
        }

        Axpy( F(-1), A, ALast );
        frobNormADiff = HermitianNorm( uplo, ALast, FROBENIUS_NORM );
    }
    while( frobNormADiff > cubeRootTol );

    MakeHermitian( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

template<typename F>
int HermitianHalley
( UpperOrLower uplo, DistMatrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("Halley");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Height must equal width");
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int height = A.Height();
    const R oneHalf = R(1)/R(2);
    const R oneThird = R(1)/R(3);

    const R epsilon = lapack::MachineEpsilon<R>();
    const R tol = 5*epsilon;
    const R cubeRootTol = Pow(tol,oneThird);
    const R a = 3;
    const R b = 1;
    const R c = 3;

    // Form the first iterate
    Scale( 1/upperBound, A );

    int numIts=0;
    R frobNormADiff;
    DistMatrix<F> ALast( g );
    DistMatrix<F> Q( 2*height, height, g );
    DistMatrix<F> QT(g), QB(g);
    PartitionDown( Q, QT,
                      QB, height );
    DistMatrix<F> C( g );
    DistMatrix<F> ATemp( g );
    do
    {
        if( numIts > 100 )
            throw std::runtime_error("Halley iteration did not converge");
        ++numIts;
        ALast = A;

        // TODO: Come up with a test for when we can use the Cholesky approach
        if( true )
        {
            //
            // The standard QR-based algorithm
            //
            MakeHermitian( uplo, A );
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            ExplicitQR( Q );
            Trrk
            ( uplo, NORMAL, ADJOINT, F(a-b/c)/Sqrt(c), QT, QB, F(b/c), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            MakeHermitian( uplo, A );
            Identity( height, height, C );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( b/c, A );
            Axpy( a-b/c, ATemp, A );
        }

        Axpy( F(-1), A, ALast );
        frobNormADiff = HermitianNorm( uplo, ALast, FROBENIUS_NORM );
    }
    while( frobNormADiff > cubeRootTol );

    MakeHermitian( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

} // namespace elem
