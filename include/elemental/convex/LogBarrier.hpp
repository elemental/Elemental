/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CONVEX_LOGBARRIER_HPP
#define CONVEX_LOGBARRIER_HPP

#include "elemental/lapack-like/Determinant.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type 
LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F>
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F> 
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

template<typename F> 
inline typename Base<F>::type
LogBarrier( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type barrier = -safeDet.kappa*safeDet.n;
#ifndef RELEASE
    PopCallStack();
#endif
    return barrier;
}

} // namespace elem

#endif // ifndef CONVEX_LOGBARRIER_HPP
