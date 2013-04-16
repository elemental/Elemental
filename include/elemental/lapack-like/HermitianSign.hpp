/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANSIGN_HPP
#define LAPACK_HERMITIANSIGN_HPP

#ifdef HAVE_PMRRR

#include "elemental/lapack-like/HermitianFunction.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"

namespace elem {

//
// Set each eigenvalue of A to its absolute value.
//

template<typename F>
inline void
HermitianSign
( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianSign");
#endif
    typedef typename Base<F>::type R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    // Compute the two-norm of A as the maximum absolute value of its eigvals
    const R twoNorm = MaxNorm( w );

    // Set the tolerance equal to n ||A||_2 eps, and invert values above it
    const int n = A.Height();
    const R eps = lapack::MachineEpsilon<R>();
    const R tolerance = n*twoNorm*eps;
    const int numLocalEigs = w.LocalHeight();
    for( int iLocal=0; iLocal<numLocalEigs; ++iLocal )
    {
        const R omega = w.GetLocal(iLocal,0);
        if( Abs(omega) < tolerance )
            w.SetLocal(iLocal,0,0);
        else if( omega > 0 )
            w.SetLocal(iLocal,0,R(1));
        else
            w.SetLocal(iLocal,0,R(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    hermitian_function::ReformHermitianMatrix( uplo, A, w, Z );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifdef HAVE_PMRRR

#endif // ifndef LAPACK_HERMITIANSIGN_HPP
