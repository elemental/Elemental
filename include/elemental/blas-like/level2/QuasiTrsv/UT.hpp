/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUASITRSV_UT_HPP
#define ELEM_QUASITRSV_UT_HPP

#include ELEM_AXPY_INC
#include ELEM_ZEROS_INC
#include ELEM_GEMV_INC

namespace elem {
namespace internal {

template<typename F>
inline void
QuasiTrsvUTUnb( Orientation orientation, const Matrix<F>& U, Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUTUnb");
        if( U.Height() != U.Width() )
            LogicError("L must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    LogicError("Not yet written");
}

template<typename F>
inline void
QuasiTrsvUT( Orientation orientation, const Matrix<F>& U, Matrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUT");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    LogicError("Not yet written");
}

template<typename F>
inline void
QuasiTrsvUT( Orientation orientation, const DistMatrix<F>& U, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::QuasiTrsvUT");
        if( U.Grid() != x.Grid() )
            LogicError("{U,x} must be distributed over the same grid");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal QuasiTrsvUT");
        if( orientation == NORMAL )
            LogicError("Invalid orientation");
    )
    LogicError("Not yet written");
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_QUASITRSV_UT_HPP
