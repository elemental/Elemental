/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRTRMM_HPP
#define BLAS_TRTRMM_HPP

#include "./Trtrmm/Unblocked.hpp"
#include "./Trtrmm/LVar1.hpp"
#include "./Trtrmm/UVar1.hpp"

namespace elem {

namespace internal {

template<typename T>
inline void
LocalTrtrmm
( Orientation orientation, UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTrtrmm");
#endif
    Trtrmm( orientation, uplo, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Trtrmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Trtrmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_TRTRMM_HPP
