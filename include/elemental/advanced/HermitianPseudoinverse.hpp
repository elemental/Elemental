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

#ifndef WITHOUT_PMRRR

namespace elemental {

namespace hermitian_pseudoinverse {

template<typename R>
class Functor
{
    const R tolerance_;
public:
    Functor( R tolerance ) : tolerance_(tolerance) { }

    R operator()( R alpha ) const
    {
        if( Abs(alpha) < tolerance_ )    
            return 0;
        else
            return 1/alpha;
    }
};

} // namespace hermitian_pseudoinverse

//
// Invert the sufficiently large eigenvalues of A.
//

template<typename R>
inline void
HermitianPseudoinverse( UpperOrLower uplo, DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianPseudoinverse");
#endif
    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<R,MC,MR> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Compute the two-norm of A as the maximum absolute value 
    // of its eigenvalues
    R maxLocalAbsEig = 0;
    const int localHeight = w.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        maxLocalAbsEig = 
            std::max(maxLocalAbsEig,Abs(w.GetLocalEntry(iLocal,0)));
    R twoNorm;
    mpi::AllReduce( &maxLocalAbsEig, &twoNorm, 1, mpi::MAX, g.VCComm() );

    // Set the tolerance equal to n ||A||_2 eps
    const int n = A.Height();
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = n*twoNorm*eps;

    // Form the pseudoinverse
    hermitian_pseudoinverse::Functor<R> f( tolerance );
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z, f );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
HermitianPseudoinverse
( UpperOrLower uplo, DistMatrix<std::complex<R>,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianPseudoinverse");
#endif
    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<std::complex<R>,MC,MR> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Compute the two-norm of A as the maximum absolute value of its
    // eigenvalues
    R maxLocalAbsEig = 0;
    const int localHeight = w.LocalHeight();
    for( int iLocal=0; iLocal<localHeight; ++iLocal )
        maxLocalAbsEig = 
            std::max(maxLocalAbsEig,Abs(w.GetLocalEntry(iLocal,0)));
    R twoNorm;
    mpi::AllReduce( &maxLocalAbsEig, &twoNorm, 1, mpi::MAX, g.VCComm() );

    // Set the tolerance equal to n ||A||_2 eps
    const int n = A.Height();
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = n*twoNorm*eps;

    // Form the pseudoinverse
    hermitian_pseudoinverse::Functor<R> f( tolerance );
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z, f );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental

#endif // WITHOUT_PMRRR
