/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRTRMM_HPP
#define ELEM_BLAS_TRTRMM_HPP

#include "./Trtrmm/Unblocked.hpp"
#include "./Trtrmm/LVar1.hpp"
#include "./Trtrmm/UVar1.hpp"

namespace elem {

template<typename T>
inline void
LocalTrtrmm
( Orientation orientation, UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalTrtrmm");
#endif
    Trtrmm( orientation, uplo, A.Matrix() );
}

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Trtrmm");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
}

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Trtrmm");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_TRTRMM_HPP
