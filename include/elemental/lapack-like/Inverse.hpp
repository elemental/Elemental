/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_INVERSE_HPP
#define ELEM_LAPACK_INVERSE_HPP

#include "elemental/blas-like/level2/ApplySymmetricPivots.hpp"
#include "elemental/lapack-like/LDL.hpp"

#include "elemental/lapack-like/Inverse/CholeskyLVar2.hpp"
#include "elemental/lapack-like/Inverse/CholeskyUVar2.hpp"
#include "elemental/lapack-like/Inverse/LUPartialPiv.hpp"

namespace elem {

template<typename F> 
inline void
Inverse( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Inverse");
#endif
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
}

template<typename F>
inline void
SymmetricInverse
( UpperOrLower uplo, Matrix<F>& A, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry cse("SymmetricInverse");
#endif
    if( uplo == LOWER )
    {
        Matrix<Int> p;
        Matrix<F> dSub;
        ldl::Pivoted( A, dSub, p, conjugate, pivotType );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );
        ApplySymmetricPivots( LOWER, A, p, conjugate );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
SymmetricInverse
( UpperOrLower uplo, DistMatrix<F>& A, bool conjugate=false,
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry cse("SymmetricInverse");
#endif
    if( uplo == LOWER )
    {
        DistMatrix<Int,VC,STAR> p( A.Grid() );
        DistMatrix<F,MD,STAR> dSub( A.Grid() );
        ldl::Pivoted( A, dSub, p, conjugate, pivotType );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );
        ApplySymmetricPivots( LOWER, A, p, conjugate );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianInverse");
#endif
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry cse("HermitianInverse");
#endif
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F> 
inline void
Inverse( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Inverse");
#endif
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
}

template<typename F>
inline void
LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalInverse");
#endif
    Inverse( A.Matrix() );
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalHPDInverse");
#endif
    HPDInverse( uplo, A.Matrix() );
}

template<typename F>
inline void
LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalSymmetricInverse");
#endif
    SymmetricInverse( uplo, A.Matrix(), conjugate, pivotType );
}

template<typename F>
inline void
LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalHermitianInverse");
#endif
    SymmetricInverse( uplo, A.Matrix(), true, pivotType );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_INVERSE_HPP
