/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_HPP
#define LAPACK_CHOLESKY_HPP

// TODO: Reorganize Cholesky implementation?
namespace elem {
template<typename F>
void LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );
} // namespace elem

#include "./Cholesky/LVar3.hpp"
#include "./Cholesky/LVar3Square.hpp"
#include "./Cholesky/UVar3.hpp"
#include "./Cholesky/UVar3Square.hpp"
#include "./Cholesky/SolveAfter.hpp"

namespace elem {

template<typename F>
inline void
LocalCholesky
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("LocalCholesky");
#endif
    Cholesky( uplo, A.Matrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::CholeskyLVar3( A );
    else
        internal::CholeskyUVar3( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
#endif
    const Grid& g = A.Grid();

    if( g.Height() == g.Width() )
    {
        if( uplo == LOWER )
            internal::CholeskyLVar3Square( A );
        else
            internal::CholeskyUVar3Square( A );
    }
    else
    {
        if( uplo == LOWER )
            internal::CholeskyLVar3( A );
        else
            internal::CholeskyUVar3( A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_CHOLESKY_HPP
