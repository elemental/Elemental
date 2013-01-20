/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_CAUCHYLIKE_HPP
#define MATRICES_CAUCHYLIKE_HPP

namespace elem {

template<typename F> 
inline void
CauchyLike
( const std::vector<F>& r, const std::vector<F>& s,
  const std::vector<F>& x, const std::vector<F>& y, 
  Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("CauchyLike");
#endif
    const int m = r.size();
    const int n = s.size();
    if( x.size() != (unsigned)m )
        throw std::logic_error("x vector was the wrong length");
    if( y.size() != (unsigned)n )
        throw std::logic_error("y vector was the wrong length");
    A.ResizeTo( m, n );

    for( int j=0; j<n; ++j )
    {
        for( int i=0; i<m; ++i )
        {
#ifndef RELEASE
            // TODO: Use tolerance instead?
            if( x[i] == y[j] )
            {
                std::ostringstream msg;
                msg << "x[" << i << "] = y[" << j << "] (" << x[i] 
                    << ") is not allowed for Cauchy-like matrices";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            A.Set( i, j, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F,Distribution U,Distribution V>
inline void
CauchyLike
( const std::vector<F>& r, const std::vector<F>& s, 
  const std::vector<F>& x, const std::vector<F>& y, 
  DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("CauchyLike");
#endif
    const int m = r.size();
    const int n = s.size();
    if( x.size() != (unsigned)m )
        throw std::logic_error("x vector was the wrong length");
    if( y.size() != (unsigned)n )
        throw std::logic_error("y vector was the wrong length");
    A.ResizeTo( m, n );

    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const int j = rowShift + jLocal*rowStride;
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const int i = colShift + iLocal*colStride;
#ifndef RELEASE
            // TODO: Use tolerance instead?
            if( x[i] == y[j] )
            {
                std::ostringstream msg;
                msg << "x[" << i << "] = y[" << j << "] (" << x[i] 
                    << ") is not allowed for Cauchy-like matrices";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            A.SetLocal( iLocal, jLocal, r[i]*s[j]/(x[i]-y[j]) );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef MATRICES_CAUCHYLIKE_HPP
