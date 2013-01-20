/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HPSDCHOLESKY_HPP
#define LAPACK_HPSDCHOLESKY_HPP

#ifndef WITHOUT_PMRRR

#include "elemental/lapack-like/HPSDSquareRoot.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/lapack-like/QR.hpp"

namespace elem {

//
// Compute the Cholesky factor of a potentially singular Hermitian semi-definite
// matrix.
//

template<typename R>
inline void
HPSDCholesky( UpperOrLower uplo, DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HPSDCholesky");
#endif
    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    if( uplo == LOWER )
    {
        LQ( A );
        MakeTrapezoidal( LEFT, LOWER, 0, A );
    }
    else
    {
        QR( A );
        MakeTrapezoidal( RIGHT, UPPER, 0, A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
HPSDCholesky( UpperOrLower uplo, DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("HPSDCholesky");
#endif
    HPSDSquareRoot( uplo, A );
    MakeHermitian( uplo, A );

    const Grid& g = A.Grid();
    if( uplo == LOWER )
    {
        DistMatrix<Complex<R>,MD,STAR> t(g);
        LQ( A, t );
        MakeTrapezoidal( LEFT, LOWER, 0, A );
    }
    else
    {
        DistMatrix<Complex<R>,MD,STAR> t(g);
        QR( A, t );
        MakeTrapezoidal( RIGHT, UPPER, 0, A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // WITHOUT_PMRRR

#endif // ifndef LAPACK_HPSDCHOLESKY_HPP
