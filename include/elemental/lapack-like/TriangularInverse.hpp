/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_TRIANGULARINVERSE_HPP
#define LAPACK_TRIANGULARINVERSE_HPP

namespace elem {
namespace internal {
template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );
} // namespace internal
} // namespace elem

#include "./TriangularInverse/LVar3.hpp"
#include "./TriangularInverse/UVar3.hpp"

namespace elem {

namespace internal {

template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTriangularInverse");
#endif
    TriangularInverse( uplo, diag, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TriangularInverseVar3
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("internal::TriangularInverseVar3");
#endif
    if( uplo == LOWER )
        TriangularInverseLVar3( diag, A );
    else
        TriangularInverseUVar3( diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TriangularInverseVar3
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("internal::TriangularInverseVar3");
#endif
    if( uplo == LOWER )
        TriangularInverseLVar3( diag, A );
    else
        TriangularInverseUVar3( diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename F>
inline void
TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
#endif
    internal::TriangularInverseVar3( uplo, diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
#endif
    internal::TriangularInverseVar3( uplo, diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_TRIANGULARINVERSE_HPP
