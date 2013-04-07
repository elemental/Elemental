/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_POLAR_HALLEY_HPP
#define LAPACk_POLAR_HALLEY_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/MakeHermitian.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/matrices/Identity.hpp"

namespace elem {

//
// A simplification of QDWH.hpp, which is loosely based on Yuji Nakatsukasa's 
// implementation of a QR-based dynamically weighted Halley iteration for the 
// polar decomposition. 
//     http://www.mathworks.com/matlabcentral/fileexchange/36830
//

namespace polar {

template<typename F>
inline int 
Halley( Matrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("polar::Halley");
#endif
    typedef typename Base<F>::type R;
    const int height = A.Height();
    const int width = A.Width();
    const R oneHalf = R(1)/R(2);
    const R oneThird = R(1)/R(3);

    if( height < width )
        throw std::logic_error("Height cannot be less than width");

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
    Matrix<F> Q( height+width, width );
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
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q );
            Gemm( NORMAL, ADJOINT, F(a-b/c)/Sqrt(c), QT, QB, F(b/c), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            Identity( width, width, C );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( b/c, A );
            Axpy( a-b/c, ATemp, A );
        }

        Axpy( F(-1), A, ALast );
        frobNormADiff = FrobeniusNorm( ALast );
    }
    while( frobNormADiff > cubeRootTol );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

template<typename F>
inline int 
Halley
( DistMatrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("polar::Halley");
#endif
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int height = A.Height();
    const int width = A.Width();
    const R oneHalf = R(1)/R(2);
    const R oneThird = R(1)/R(3);

    if( height < width )
        throw std::logic_error("Height cannot be less than width");

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
    DistMatrix<F> Q( height+width, width, g );
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
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q );
            Gemm( NORMAL, ADJOINT, F(a-b/c)/Sqrt(c), QT, QB, F(b/c), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            Identity( width, width, C );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( b/c, A );
            Axpy( a-b/c, ATemp, A );
        }

        Axpy( F(-1), A, ALast );
        frobNormADiff = FrobeniusNorm( ALast );
    }
    while( frobNormADiff > cubeRootTol );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

} // namespace polar

namespace hermitian_polar {

template<typename F>
inline int
Halley( UpperOrLower uplo, Matrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("hermitian_polar::Halley");
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
            qr::Explicit( Q );
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
        frobNormADiff = HermitianFrobeniusNorm( uplo, ALast );
    }
    while( frobNormADiff > cubeRootTol );

    MakeHermitian( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

template<typename F>
inline int
Halley( UpperOrLower uplo, DistMatrix<F>& A, typename Base<F>::type upperBound )
{
#ifndef RELEASE
    PushCallStack("hermitian_polar::Halley");
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
            qr::Explicit( Q );
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
        frobNormADiff = HermitianFrobeniusNorm( uplo, ALast );
    }
    while( frobNormADiff > cubeRootTol );

    MakeHermitian( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return numIts;
}

} // namespace hermitian_polar

} // namespace elem

#endif // ifndef LAPACK_POLAR_HALLEY_HPP
