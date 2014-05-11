/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CAUCHY_HPP
#define ELEM_CAUCHY_HPP

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

template<typename F1,typename F2,Dist U,Dist V>
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
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
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

template<typename F1,typename F2,Dist U,Dist V>
inline void
Cauchy
( BlockDistMatrix<F1,U,V>& A, const std::vector<F2>& x, const std::vector<F2>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Cauchy"))
    const Int m = x.size();
    const Int n = y.size();
    A.Resize( m, n );

    const F1 one = F1(1);
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
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

} // namespace elem

#endif // ifndef ELEM_CAUCHY_HPP
