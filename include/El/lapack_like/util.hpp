/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UTIL_HPP
#define EL_UTIL_HPP

namespace El {

// Median
// ======
template<typename Real>
ValueInt<Real> Median( const Matrix<Real>& x );
template<typename Real>
ValueInt<Real> Median( const AbstractDistMatrix<Real>& x );

// Sort
// ====
template<typename Real>
void Sort( Matrix<Real>& X, SortType sort=ASCENDING );
template<typename Real>
void Sort( AbstractDistMatrix<Real>& X, SortType sort=ASCENDING );

template<typename Real>
vector<ValueInt<Real>> TaggedSort
( const Matrix<Real>& x, SortType sort=ASCENDING );
template<typename Real>
vector<ValueInt<Real>> TaggedSort
( const AbstractDistMatrix<Real>& x, SortType sort=ASCENDING );

} // namespace El

#endif // ifndef EL_UTIL_HPP
