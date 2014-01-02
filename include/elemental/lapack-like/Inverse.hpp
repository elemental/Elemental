/*
   Copyright (c) 2009-2014, Jack Poulson
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
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDInverse"))
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
    DEBUG_ONLY(CallStackEntry cse("SymmetricInverse"))
    if( uplo == LOWER )
    {
        Matrix<Int> p;
        Matrix<F> dSub;
        ldl::Pivoted( A, dSub, p, conjugate, pivotType );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );
        ApplyInverseSymmetricPivots( LOWER, A, p, conjugate );
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
    DEBUG_ONLY(CallStackEntry cse("SymmetricInverse"))
    if( uplo == LOWER )
    {
        DistMatrix<Int,VC,STAR> p( A.Grid() );
        DistMatrix<F,MD,STAR> dSub( A.Grid() );
        ldl::Pivoted( A, dSub, p, conjugate, pivotType );
        TriangularInverse( LOWER, UNIT, A ); 
        Trdtrmm( LOWER, A, dSub, conjugate );
        ApplyInverseSymmetricPivots( LOWER, A, p, conjugate );
    }
    else
        LogicError("This option is not yet supported");
}

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
inline void
HermitianInverse
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F> 
inline void
Inverse( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDInverse"))
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
}

template<typename F>
inline void
LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalInverse"))
    Inverse( A.Matrix() );
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHPDInverse"))
    HPDInverse( uplo, A.Matrix() );
}

template<typename F>
inline void
LocalSymmetricInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, bool conjugate=false, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalSymmetricInverse"))
    SymmetricInverse( uplo, A.Matrix(), conjugate, pivotType );
}

template<typename F>
inline void
LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, 
  LDLPivotType pivotType=BUNCH_KAUFMAN_A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHermitianInverse"))
    SymmetricInverse( uplo, A.Matrix(), true, pivotType );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_INVERSE_HPP
