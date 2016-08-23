/* ========================================================================== */
/* === ldl.hpp: include file for the LDL package ============================ */
/* ========================================================================== */

/* 
 * Copyright (c) Timothy A Davis, http://www.suitesparse.com.
 * All Rights Reserved. 
 *
 * Copyright (c) Jack Poulson, https://github.com/elemental/Elemental.
 * All Rights Reserved.
 *
 * Your use or distribution of LDL or any modified version of
 * LDL implies that you agree to this License.

 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
 * USA
 *
 * Permission is hereby granted to use or copy this program under the
 * terms of the GNU LGPL, provided that the Copyright, this License,
 * and the Availability of the original version is retained on all copies.
 * User documentation of any code that uses this code or any modified
 * version of this code must cite the Copyright, this License, the
 * Availability note, and "Used by permission." Permission to modify
 * the code and to distribute modified code is granted, provided the
 * Copyright, this License, and the Availability note are retained,
 * and a notice that the code was modified is included.
 */
#ifndef EL_SUITE_SPARSE_LDL_HPP
#define EL_SUITE_SPARSE_LDL_HPP

#include "ElSuiteSparse/config.h"
#include <complex>

namespace suite_sparse {
namespace ldl {

using std::complex;

template<typename IntType=int>
void Symbolic
( IntType n,
  const IntType* Ap,
  const IntType* Ai,
        IntType* Lp,
        IntType* Parent,
        IntType* Lnz,
        IntType* Flag, 
  const IntType* P,
        IntType* Pinv );

template<typename F=double,typename IntType=int>
IntType Numeric
( IntType n, 
  const IntType* Ap,
  const IntType* Ai,
  const F* Ax,
  const IntType* Lp,
  const IntType* Parent,
        IntType* Lnz,
        IntType* Li,
        F* Lx, 
        F* D,
        F* Y,
        IntType* Pattern,
        IntType* Flag,
  const IntType* P,
  const IntType* Pinv,
  bool conjugate=false );

template<typename F=double,typename IntType=int>
void LSolve
( IntType m, 
        F* X,
  const IntType* Lp,
  const IntType* Li,
  const F* Lx ) ;
template<typename F=double,typename IntType=int>
void LSolveMulti
( bool onLeft,
  IntType m, 
  IntType n,
        F* X,
  IntType XLDim,
  const IntType* Lp,
  const IntType* Li,
  const F* Lx ) ;

template<typename F=double,typename IntType=int>
void DSolve( IntType n, F* X, const F* D );
template<typename F=double,typename IntType=int>
void DSolveMulti
( bool onLeft, IntType m, IntType n, F* X, IntType XLDim, const F* D );

template<typename F=double,typename IntType=int>
void LTSolve
( IntType m,
        F* X,
  const IntType* Lp,
  const IntType* Li,
  const F* Lx,
  bool conjugate=false );
template<typename F=double,typename IntType=int>
void LTSolveMulti
( bool onLeft,
  IntType m,
  IntType n,
        F* X,
  IntType XLDim,
  const IntType* Lp,
  const IntType* Li,
  const F* Lx,
  bool conjugate=false );

template<typename F=double,typename IntType=int>
void Perm( IntType n, F* X, const F* B, const IntType* P );

template<typename F=double,typename IntType=int>
void PermT( IntType n, F* X, const F* B, const IntType* P );

// TODO: Multi versions of Perm and PermT?

template<typename IntType=int>
IntType ValidPerm( IntType n, const IntType* P, IntType* Flag );

template<typename IntType=int>
IntType ValidMatrix( IntType n, const IntType* Ap, const IntType* Ai );

} // namespace ldl
} // namespace suite_sparse

#define EL_LDL_DATE "Dec 27, 2015"
#define EL_LDL_VERSION_CODE(main,sub) ((main) * 1000 + (sub))
#define EL_LDL_MAIN_VERSION 2
#define EL_LDL_SUB_VERSION 2
#define EL_LDL_SUBSUB_VERSION 1
#define EL_LDL_VERSION EL_LDL_VERSION_CODE( \
 EL_LDL_MAIN_VERSION,EL_LDL_SUB_VERSION)

#include "./ldl/impl.hpp"

#endif // ifndef EL_SUITE_SPARSE_LDL_HPP
