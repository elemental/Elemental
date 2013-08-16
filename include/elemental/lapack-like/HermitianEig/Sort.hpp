/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_HERMITIANEIG_SORT_HPP
#define ELEM_LAPACK_HERMITIANEIG_SORT_HPP

namespace elem {

namespace hermitian_eig {

template<typename R>
inline void
Sort( Matrix<R>& w, bool ascending=true )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_eig::Sort");
#endif
    R* wBuffer = w.Buffer();
    const Int k = w.Height();
    if( ascending )
        std::sort( &wBuffer[0], &wBuffer[k] );
    else
        std::sort( &wBuffer[0], &wBuffer[k], std::greater<R>() );
}

template<typename R>
inline void
Sort( DistMatrix<R,VR,STAR>& w, bool ascending=true )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_eig::Sort");
#endif
    const Int k = w.Height();

    // Gather a full copy of w on each process and locally sort
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );
    R* wBuffer = w_STAR_STAR.Buffer();
    if( ascending )
        std::sort( &wBuffer[0], &wBuffer[k] );
    else
        std::sort( &wBuffer[0], &wBuffer[k], std::greater<R>() );

    // Refill the distributed w with the sorted values
    w = w_STAR_STAR;
}

template<typename F>
inline void
Sort( Matrix<BASE(F)>& w, Matrix<F>& Z, bool ascending=true )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_eig::Sort");
#endif
    typedef BASE(F) R;
    const Int n = Z.Height();
    const Int k = Z.Width();

    // Initialize the pairs of indices and eigenvalues
    std::vector<ValueInt<R> > pairs( k );
    for( Int i=0; i<k; ++i )
    {
        pairs[i].value = w.Get(i,0);
        pairs[i].index = i;
    }

    // Sort the eigenvalues and simultaneously form the permutation
    if( ascending )
        std::sort( pairs.begin(), pairs.end(), ValueInt<R>::Lesser );
    else
        std::sort( pairs.begin(), pairs.end(), ValueInt<R>::Greater );

    // Reorder the eigenvectors and eigenvalues using the new ordering
    Matrix<F> ZPerm( n, k );
    for( Int j=0; j<k; ++j )
    {
        const Int source = pairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), n );
        w.Set(j,0,pairs[j].value);
    }
    Z = ZPerm;
}

template<typename F>
inline void
Sort( DistMatrix<BASE(F),VR,STAR>& w, DistMatrix<F>& Z, bool ascending=true )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_eig::Sort");
#endif
    typedef BASE(F) R;
    const Int n = Z.Height();
    const Int k = Z.Width();
    const Grid& g = Z.Grid();

    DistMatrix<F,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );

    // Initialize the pairs of indices and eigenvalues
    std::vector<ValueInt<R> > pairs( k );
    for( Int i=0; i<k; ++i )
    {
        pairs[i].value = w_STAR_STAR.GetLocal(i,0);
        pairs[i].index = i;
    }

    // Sort the eigenvalues and simultaneously form the permutation
    if( ascending )
        std::sort( pairs.begin(), pairs.end(), ValueInt<R>::Lesser );
    else
        std::sort( pairs.begin(), pairs.end(), ValueInt<R>::Greater );

    // Locally reorder the eigenvectors and eigenvalues using the new ordering
    const Int nLocal = Z_VC_STAR.LocalHeight();
    DistMatrix<F,VC,STAR> ZPerm_VC_STAR( n, k, g );
    for( Int j=0; j<k; ++j )
    {
        const Int source = pairs[j].index;
        MemCopy
        ( ZPerm_VC_STAR.Buffer(0,j), 
          Z_VC_STAR.LockedBuffer(0,source), nLocal );
        w_STAR_STAR.SetLocal(j,0,pairs[j].value);
    }
    Z_VC_STAR.Empty();

    Z = ZPerm_VC_STAR;
    w = w_STAR_STAR;
}

} // namespace hermitian_eig

} // namespace elem

#endif // ifndef ELEM_LAPACK_HERMITIANEIG_SORT_HPP
