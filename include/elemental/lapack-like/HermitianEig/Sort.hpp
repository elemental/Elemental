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

#include "elemental/lapack-like/Sort.hpp"

namespace elem {

namespace hermitian_eig {

template<typename F>
inline void
Sort( Matrix<BASE(F)>& w, Matrix<F>& Z, SortType sort=ASCENDING )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::Sort");
#endif
    if( sort == UNSORTED )
        return;

    // Initialize the pairs of indices and eigenvalues
    typedef Base<F> Real;
    std::vector<ValueInt<Real>> pairs = TaggedSort( w, sort );

    // Reorder the eigenvectors and eigenvalues using the new ordering
    const Int n = Z.Height();
    const Int k = Z.Width();
    Matrix<F> ZPerm( n, k );
    for( Int j=0; j<k; ++j )
    {
        const Int source = pairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), n );
        w.Set(j,0,pairs[j].value);
    }
    Z = ZPerm;
}

template<typename F,Distribution U,Distribution V>
inline void
Sort( DistMatrix<BASE(F),U,V>& w, DistMatrix<F>& Z, SortType sort=ASCENDING )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::Sort");
#endif
    if( sort == UNSORTED )
        return;

    // Get the sorted eigenvalue information
    typedef Base<F> Real;
    std::vector<ValueInt<Real>> pairs = TaggedSort( w, sort );

    // Locally reorder the eigenvectors and eigenvalues using the new ordering
    const Int n = Z.Height();
    const Int k = Z.Width();
    const Grid& g = Z.Grid();
    DistMatrix<F,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<F,VC,STAR> ZPerm_VC_STAR( n, k, g );
    const Int nLocal = Z_VC_STAR.LocalHeight();
    for( Int j=0; j<k; ++j )
    {
        MemCopy
        ( ZPerm_VC_STAR.Buffer(0,j), 
          Z_VC_STAR.LockedBuffer(0,pairs[j].index), nLocal );
        w.Set( j, 0, pairs[j].value );
    }
    Z_VC_STAR.Empty();
    Z = ZPerm_VC_STAR;
}

} // namespace hermitian_eig

} // namespace elem

#endif // ifndef ELEM_LAPACK_HERMITIANEIG_SORT_HPP
