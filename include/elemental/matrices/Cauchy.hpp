/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_CAUCHY_HPP
#define MATRICES_CAUCHY_HPP

namespace elem {

template<typename F> 
inline void
Cauchy( Matrix<F>& A, const std::vector<F>& x, const std::vector<F>& y )
{
#ifndef RELEASE
    CallStackEntry entry("Cauchy");
#endif
    const int m = x.size();
    const int n = y.size();
    A.ResizeTo( m, n );

    const F one = F(1);
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
                    << ") is not allowed for Cauchy matrices";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            A.Set( i, j, one/(x[i]-y[j]) );
        }
    }
}

template<typename F,Distribution U,Distribution V>
inline void
Cauchy( DistMatrix<F,U,V>& A, const std::vector<F>& x, const std::vector<F>& y )
{
#ifndef RELEASE
    CallStackEntry entry("Cauchy");
#endif
    const int m = x.size();
    const int n = y.size();
    A.ResizeTo( m, n );

    const F one = F(1);
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    const int colShift = A.ColShift();
    const int rowShift = A.RowShift();
    const int colStride = A.ColStride();
    const int rowStride = A.RowStride();
    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*rowStride;
        for( int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const int i = colShift + iLoc*colStride;
#ifndef RELEASE
            // TODO: Use tolerance instead?
            if( x[i] == y[j] )
            {
                std::ostringstream msg;
                msg << "x[" << i << "] = y[" << j << "] (" << x[i] 
                    << ") is not allowed for Cauchy matrices";
                throw std::logic_error( msg.str().c_str() );
            }
#endif
            A.SetLocal( iLoc, jLoc, one/(x[i]-y[j]) );
        }
    }
}

} // namespace elem

#endif // ifndef MATRICES_CAUCHY_HPP
