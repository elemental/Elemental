/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HERMITIANTRIDIAGEIG_SORT_HPP
#define ELEM_HERMITIANTRIDIAGEIG_SORT_HPP

#include ELEM_SORT_INC

namespace elem {

// This routine is slightly more general than necessary (complex support) so
// that it may also be used for sorting Hermitian eigenpairs
namespace herm_eig {

template<typename F>
inline void
Sort( Matrix<Base<F>>& w, Matrix<F>& Z, SortType sort=ASCENDING )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::Sort"))
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

template<typename F,Dist U1,Dist V1,
                    Dist U2,Dist V2>
inline void
Sort
( DistMatrix<Base<F>,U1,V1>& w, DistMatrix<F,U2,V2>& Z, 
  SortType sort=ASCENDING )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::Sort"))
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

} // namespace herm_eig

} // namespace elem

#endif // ifndef ELEM_HERMITIANTRIDIAGEIG_SORT_HPP
