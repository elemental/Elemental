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
// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.
//

template<typename F>
inline void
Polar( Matrix<F>& A, Matrix<F>& P )
{
#ifndef RELEASE
    PushCallStack("Polar");
#endif
    typedef typename Base<F>::type R;
    const int n = A.Width();

    // Get the SVD of A
    Matrix<R> s;
    Matrix<F> U, V;
    U = A;
    SVD( U, s, V );

    // Form Q := U V^H in A
    MakeZeros( A );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    // Form P := V Sigma V^H in P
    Zeros( n, n, P );
    hermitian_function::ReformHermitianMatrix( LOWER, P, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Polar( DistMatrix<F>& A, DistMatrix<F>& P )
{
#ifndef RELEASE
    PushCallStack("Polar");
#endif
    typedef typename Base<F>::type R;
    const Grid& g = A.Grid();
    const int n = A.Width();

    // Get the SVD of A
    DistMatrix<R,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    SVD( U, s, V );

    // Form Q := U V^H in A
    MakeZeros( A );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    // Form P := V Sigma V^H in P
    Zeros( n, n, P );
    hermitian_function::ReformHermitianMatrix( LOWER, P, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
