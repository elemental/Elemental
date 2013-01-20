/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_TRACE_HPP
#define LAPACK_TRACE_HPP

namespace elem {

template<typename F> 
inline F Trace( const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Trace");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot compute trace of nonsquare matrix");
    const Grid& g = A.Grid();

    DistMatrix<F,MD,STAR> d(g);
    A.GetDiagonal( d );
    F localTrace = 0;
    if( d.Participating() )
    {
        const int nLocalDiag = d.LocalHeight();
        for( int iLocal=0; iLocal<nLocalDiag; ++iLocal )
            localTrace += d.GetLocal(iLocal,0);
    }
    F trace;
    mpi::AllReduce( &localTrace, &trace, 1, mpi::SUM, g.VCComm() );
#ifndef RELEASE
    PopCallStack();
#endif
    return trace;
}

template<typename F>
inline F Trace( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Trace");
#endif
    if( A.Height() != A.Width() )
        throw std::logic_error("Cannot compute trace of nonsquare matrix");

    Matrix<F> d;
    A.GetDiagonal( d );
    F trace = 0;
    const int n = A.Height();
    for( int i=0; i<n; ++i )
        trace += d.Get(i,0);
#ifndef RELEASE
    PopCallStack();
#endif
    return trace;
}

} // namespace elem

#endif // ifndef LAPACK_TRACE_HPP
