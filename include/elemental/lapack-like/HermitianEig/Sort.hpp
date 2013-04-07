/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANEIG_SORT_HPP
#define LAPACK_HERMITIANEIG_SORT_HPP

namespace elem {

template<typename R>
struct IndexValuePair {
    int index;
    R value;    

    static bool Lesser
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return a.value < b.value; }
    static bool Greater
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return a.value > b.value; }
};

template<typename R>
struct IndexValuePair<Complex<R> > {
    int index;
    Complex<R> value;    

    static bool Lesser
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return Abs(a.value) < Abs(b.value); }
    static bool Greater
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return Abs(a.value) > Abs(b.value); }
};

namespace hermitian_eig {

template<typename R>
inline void
Sort( Matrix<R>& w, bool ascending=true )
{
#ifndef RELEASE
    PushCallStack("hermitian_eig::Sort");
#endif
    R* wBuffer = w.Buffer();
    const int k = w.Height();
    if( ascending )
        std::sort( &wBuffer[0], &wBuffer[k] );
    else
        std::sort( &wBuffer[0], &wBuffer[k], std::greater<R>() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
Sort( DistMatrix<R,VR,STAR>& w, bool ascending=true )
{
#ifndef RELEASE
    PushCallStack("hermitian_eig::Sort");
#endif
    const int k = w.Height();

    // Gather a full copy of w on each process and locally sort
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );
    R* wBuffer = w_STAR_STAR.Buffer();
    if( ascending )
        std::sort( &wBuffer[0], &wBuffer[k] );
    else
        std::sort( &wBuffer[0], &wBuffer[k], std::greater<R>() );

    // Refill the distributed w with the sorted values
    w = w_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Sort( Matrix<typename Base<F>::type>& w, Matrix<F>& Z, bool ascending=true )
{
#ifndef RELEASE
    PushCallStack("hermitian_eig::Sort");
#endif
    typedef typename Base<F>::type R;
    const int n = Z.Height();
    const int k = Z.Width();

    // Initialize the pairs of indices and eigenvalues
    std::vector<IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w.Get(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    if( ascending )
        std::sort( pairs.begin(), pairs.end(), IndexValuePair<R>::Lesser );
    else
        std::sort( pairs.begin(), pairs.end(), IndexValuePair<R>::Greater );

    // Reorder the eigenvectors and eigenvalues using the new ordering
    Matrix<F> ZPerm( n, k );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), n );
        w.Set(j,0,pairs[j].value);
    }
    Z = ZPerm;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Sort
( DistMatrix<typename Base<F>::type,VR,STAR>& w, 
  DistMatrix<F>& Z, bool ascending=true )
{
#ifndef RELEASE
    PushCallStack("hermitian_eig::Sort");
#endif
    typedef typename Base<F>::type R;
    const int n = Z.Height();
    const int k = Z.Width();
    const Grid& g = Z.Grid();

    DistMatrix<F,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );

    // Initialize the pairs of indices and eigenvalues
    std::vector<IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w_STAR_STAR.GetLocal(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    if( ascending )
        std::sort( pairs.begin(), pairs.end(), IndexValuePair<R>::Lesser );
    else
        std::sort( pairs.begin(), pairs.end(), IndexValuePair<R>::Greater );

    // Locally reorder the eigenvectors and eigenvalues using the new ordering
    const int nLocal = Z_VC_STAR.LocalHeight();
    DistMatrix<F,VC,STAR> ZPerm_VC_STAR( n, k, g );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy
        ( ZPerm_VC_STAR.Buffer(0,j), 
          Z_VC_STAR.LockedBuffer(0,source), nLocal );
        w_STAR_STAR.SetLocal(j,0,pairs[j].value);
    }
    Z_VC_STAR.Empty();

    Z = ZPerm_VC_STAR;
    w = w_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace hermitian_eig

} // namespace elem

#endif // ifndef LAPACK_HERMITIANEIG_SORT_HPP
