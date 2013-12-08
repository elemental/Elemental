/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CONVEX_LOGBARRIER_HPP
#define ELEM_CONVEX_LOGBARRIER_HPP

#include "elemental/lapack-like/Determinant.hpp"

namespace elem {

template<typename F>
inline BASE(F) 
LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

#ifndef SWIG
template<typename F>
inline BASE(F)
LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry cse("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}
#endif

template<typename F> 
inline BASE(F)
LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

#ifndef SWIG
template<typename F> 
inline BASE(F)
LogBarrier( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry cse("LogBarrier");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}
#endif

} // namespace elem

#endif // ifndef ELEM_CONVEX_LOGBARRIER_HPP
