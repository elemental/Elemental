/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_CAUCHYLIKE_HPP
#define ELEM_MATRICES_CAUCHYLIKE_HPP

namespace elem {

template<typename F1,typename F2> 
inline void
CauchyLike
( Matrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s,
  const std::vector<F2>& x, const std::vector<F2>& y )
{
#ifndef RELEASE
    CallStackEntry cse("CauchyLike");
#endif
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.ResizeTo( m, n );

    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
#ifndef RELEASE
            // TODO: Use tolerance instead?
            if( x[i] == y[j] )
            {
                std::ostringstream msg;
                msg << "x[" << i << "] = y[" << j << "] (" << x[i] 
                    << ") is not allowed for Cauchy-like matrices";
                LogicError( msg.str() );
            }
#endif
            A.Set( i, j, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
}

template<typename F> 
inline Matrix<F>
CauchyLike
( const std::vector<F>& r, const std::vector<F>& s,
  const std::vector<F>& x, const std::vector<F>& y )
{
    Matrix<F> A;
    CauchyLike( A, r, s, x, y );
    return A;
}

template<typename F1,typename F2,Distribution U,Distribution V>
inline void
CauchyLike
( DistMatrix<F1,U,V>& A,
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y )
{
#ifndef RELEASE
    CallStackEntry cse("CauchyLike");
#endif
    const Int m = r.size();
    const Int n = s.size();
    if( x.size() != (Unsigned)m )
        LogicError("x vector was the wrong length");
    if( y.size() != (Unsigned)n )
        LogicError("y vector was the wrong length");
    A.ResizeTo( m, n );

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
#ifndef RELEASE
            // TODO: Use tolerance instead?
            if( x[i] == y[j] )
            {
                std::ostringstream msg;
                msg << "x[" << i << "] = y[" << j << "] (" << x[i] 
                    << ") is not allowed for Cauchy-like matrices";
                LogicError( msg.str() );
            }
#endif
            A.SetLocal( iLoc, jLoc, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
CauchyLike
( const Grid& g,
  const std::vector<F>& r, const std::vector<F>& s, 
  const std::vector<F>& x, const std::vector<F>& y )
{
    DistMatrix<F,U,V> A(g);
    CauchyLike( A, r, s, x, y );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_CAUCHYLIKE_HPP
