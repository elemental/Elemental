/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRDTRMM_HPP
#define ELEM_BLAS_TRDTRMM_HPP

namespace elem {

template<typename T>
inline void
LocalTrdtrmm
( Orientation orientation, UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalTrdtrmm");
#endif
    Trdtrmm( orientation, uplo, A.Matrix() );
}

} // namespace elem

#include "./Trdtrmm/Unblocked.hpp"
#include "./Trdtrmm/LVar1.hpp"
#include "./Trdtrmm/UVar1.hpp"

namespace elem {

template<typename F>
inline void
Trdtrmm( Orientation orientation, UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Trdtrdmm");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( orientation, A );
    else
        internal::TrdtrmmUVar1( orientation, A );
}

template<typename F>
inline void
Trdtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Trdtrmm");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( orientation, A );
    else
        internal::TrdtrmmUVar1( orientation, A );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_TRDTRMM_HPP
