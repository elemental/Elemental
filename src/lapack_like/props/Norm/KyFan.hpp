/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_KYFAN_HPP
#define EL_NORM_KYFAN_HPP

namespace El {

template<typename F> 
Base<F> KyFanNorm( const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("KyFanNorm"))
    return KyFanSchattenNorm( A, k, Base<F>(1) );
}

template<typename F>
Base<F> HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianKyFanNorm"))
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricKyFanNorm"))
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F> 
Base<F> KyFanNorm( const AbstractDistMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("KyFanNorm"))
    return KyFanSchattenNorm( A, k, Base<F>(1) );
}

template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianKyFanNorm"))
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricKyFanNorm"))
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

} // namespace El

#endif // ifndef EL_NORM_KYFAN_HPP
