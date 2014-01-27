/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRDTRMM_HPP
#define ELEM_TRDTRMM_HPP

namespace elem {

template<typename T>
inline void
LocalTrdtrmm
( UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrdtrmm"))
    Trdtrmm( uplo, A.Matrix(), conjugate );
}

template<typename T>
inline void
LocalTrdtrmm
( UpperOrLower uplo, 
  DistMatrix<T,STAR,STAR>& A, const DistMatrix<T,STAR,STAR>& dOff, 
  bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("LocalTrdtrmm"))
    Trdtrmm( uplo, A.Matrix(), dOff.LockedMatrix(), conjugate );
}

} // namespace elem

#include "./Trdtrmm/Unblocked.hpp"
#include "./Trdtrmm/LVar1.hpp"
#include "./Trdtrmm/UVar1.hpp"

namespace elem {

template<typename F>
inline void
Trdtrmm( UpperOrLower uplo, Matrix<F>& A, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrdmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( A, conjugate );
    else
        internal::TrdtrmmUVar1( A, conjugate );
}

template<typename F>
inline void
Trdtrmm
( UpperOrLower uplo, Matrix<F>& A, const Matrix<F>& dOff, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrdmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

template<typename F>
inline void
Trdtrmm( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( A, conjugate );
    else
        internal::TrdtrmmUVar1( A, conjugate );
}

template<typename F>
inline void
Trdtrmm
( UpperOrLower uplo, 
  DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dOff, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trdtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( A, dOff, conjugate );
    else
        LogicError("Not yet written");
}

} // namespace elem

#endif // ifndef ELEM_TRDTRMM_HPP
