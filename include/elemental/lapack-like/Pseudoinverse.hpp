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
// Replace A with its pseudoinverse
//

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
    R maxLocalVal= 0;
    const int numLocalVals = s.LocalHeight();
    for( int iLocal=0; iLocal<numLocalVals; ++iLocal )
    {
        const R sigma = s.GetLocal(iLocal,0);
        maxLocalVal = std::max(maxLocalVal,sigma);
    }
    R twoNorm;
    mpi::AllReduce( &maxLocalVal, &twoNorm, 1, mpi::MAX, g.VCComm() );

    // Set the tolerance equal to k ||A||_2 eps and invert above tolerance
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = k*twoNorm*eps;
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
    Gemm( NORMAL, ADJOINT, (F)1, V, U, (F)0, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
