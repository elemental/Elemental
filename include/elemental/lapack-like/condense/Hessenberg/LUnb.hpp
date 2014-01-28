/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HESSENBERG_LUNB_HPP
#define ELEM_HESSENBERG_LUNB_HPP

#include ELEM_GEMV_INC
#include ELEM_GER_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace hessenberg {

template<typename F>
inline void LUnb( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::LUnb"))
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.Resize( tHeight, 1 );

    for( Int k=0; k<n-1; ++k )
    {
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+1 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   );

        // TODO
        LogicError("This routine not yet written");
    }
}

template<typename F> 
inline void LUnb( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(CallStackEntry cse("hessenberg::LUnb"))
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Int tHeight = Max(n-1,0);
    t.SetGrid( g );
    t.Resize( tHeight, 1 );

    for( Int k=0; k<n-1; ++k )
    {
        auto a12      = ViewRange( A, k,   k+1, k+1, n   );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+1 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, n   );
        auto A22      = ViewRange( A, k+1, k+1, n,   n   );
        auto A2       = ViewRange( A, k+1, 0,   n,   n   );

        // TODO
        LogicError("This routine not yet written");
    }
}

} // namespace hessenberg
} // namespace elem

#endif // ifndef ELEM_HESSENBERG_LUNB_HPP
