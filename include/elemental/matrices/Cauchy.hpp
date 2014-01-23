/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_CAUCHY_HPP
#define ELEM_MATRICES_CAUCHY_HPP

namespace elem {

template<typename F1,typename F2> 
inline void
Cauchy( Matrix<F1>& A, const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );

    const F1 one = F1(1);
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            ) 
            A.Set( i, j, one/(x[i]-y[j]) );
        }
    }
}

template<typename F> 
inline Matrix<F>
Cauchy( const std::vector<F>& x, const std::vector<F>& y )
{
    Matrix<F> A;
    Cauchy( A, x, y );
    return A;
}

template<typename F1,typename F2,Distribution U,Distribution V>
inline void
Cauchy
( DistMatrix<F1,U,V>& A, const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );

    const F1 one = F1(1);
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            DEBUG_ONLY(
                // TODO: Use tolerance instead?
                if( x[i] == y[j] )
                    LogicError
                    ( "x[", i, "] = y[", j, "] (", x[i], 
                      ") is not allowed for Cauchy matrices" );
            )
            A.SetLocal( iLoc, jLoc, one/(x[i]-y[j]) );
        }
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Cauchy( const Grid& g, const std::vector<F>& x, const std::vector<F>& y )
{
    DistMatrix<F,U,V> A(g);
    Cauchy( A, x, y );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_CAUCHY_HPP
