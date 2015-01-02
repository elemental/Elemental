/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_SCHATTEN_HPP
#define EL_NORM_SCHATTEN_HPP

namespace El {

template<typename F> 
Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SchattenNorm"))
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSchattenNorm"))
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSchattenNorm"))
    const Int minDim = A.Height();
    return SymmetricKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F> 
Base<F> SchattenNorm( const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SchattenNorm"))
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSchattenNorm"))
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricSchattenNorm"))
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

} // namespace El

#endif // ifndef EL_NORM_SCHATTEN_HPP
