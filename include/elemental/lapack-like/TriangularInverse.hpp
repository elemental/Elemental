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
template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );
} // namespace elem

#include "elemental/lapack-like/TriangularInverse/LVar3.hpp"
#include "elemental/lapack-like/TriangularInverse/UVar3.hpp"

namespace elem {

namespace triangular_inverse {

template<typename F>
inline void
Var3( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::Var3");
#endif
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Var3( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("triangular_inverse::Var3");
#endif
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace triangular_inverse

template<typename F>
inline void
TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
#endif
    triangular_inverse::Var3( uplo, diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TriangularInverse( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
#endif
    triangular_inverse::Var3( uplo, diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("LocalTriangularInverse");
#endif
    TriangularInverse( uplo, diag, A.Matrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_TRIANGULARINVERSE_HPP
