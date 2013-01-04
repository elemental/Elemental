/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

//
// Replace A with its pseudoinverse
//

template<typename F>
inline void
Pseudoinverse( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Pseudoinverse");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::max(m,n);

    // Get the SVD of A
    Matrix<R> s;
    Matrix<F> U, V;
    U = A;
    SVD( U, s, V );

    // Compute the two-norm of A as the maximum singular value
    const R twoNorm = Norm( s, INFINITY_NORM );

    // Set the tolerance equal to k ||A||_2 eps and invert above tolerance
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = k*twoNorm*eps;
    const int numVals = s.Height();
    for( int i=0; i<numVals; ++i )
    {
        const R sigma = s.Get(i,0);
        if( sigma < tolerance )
            s.Set(i,0,0);
        else
            s.Set(i,0,1/sigma);
    }

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Zeros( n, m, A );
    Gemm( NORMAL, ADJOINT, F(1), V, U, F(0), A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Pseudoinverse( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Pseudoinverse");
#endif
    typedef typename Base<F>::type R;

    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::max(m,n);

    // Get the SVD of A
    DistMatrix<R,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    SVD( U, s, V );

    // Compute the two-norm of A as the maximum singular value
    const R twoNorm = Norm( s, INFINITY_NORM );

    // Set the tolerance equal to k ||A||_2 eps and invert above tolerance
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = k*twoNorm*eps;
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocal(iLocal,0);
        if( sigma < tolerance )
            s.SetLocal(iLocal,0,0);
        else
            s.SetLocal(iLocal,0,1/sigma);
    }

    // Scale U with the singular values, U := U Sigma
    DiagonalScale( RIGHT, NORMAL, s, U );

    // Form pinvA = (U Sigma V^H)^H = V (U Sigma)^H
    Zeros( n, m, A );
    Gemm( NORMAL, ADJOINT, F(1), V, U, F(0), A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
