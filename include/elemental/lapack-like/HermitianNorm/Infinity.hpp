/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_INFINITY_HPP
#define LAPACK_HERMITIANNORM_INFINITY_HPP

namespace elem {
namespace internal {

// The operator L1 and Linf norms for Hermitian matrices are identical. The 
// former is the maximum column L1 norm and the latter is the maximum row L1 
// norm. Hermiticity implies their equivalence.

template<typename F>
inline typename Base<F>::type
HermitianInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianInfinityNorm");
#endif
    typedef typename Base<F>::type R;
    R maxRowSum = HermitianNorm( uplo, A, ONE_NORM );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

template<typename F>
inline typename Base<F>::type
HermitianInfinityNorm
( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("internal::HermitianInfinityNorm");
#endif
    typedef typename Base<F>::type R;
    R maxRowSum = HermitianNorm( uplo, A, ONE_NORM );
#ifndef RELEASE
    PopCallStack();
#endif
    return maxRowSum;
}

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_INFINITY_HPP
