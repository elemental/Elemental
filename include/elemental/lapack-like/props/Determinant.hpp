/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DETERMINANT_HPP
#define ELEM_DETERMINANT_HPP

#include "./Determinant/Cholesky.hpp"
#include "./Determinant/LUPartialPiv.hpp"

namespace elem {

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    Matrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    DistMatrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<Base<F>> 
SafeHPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    Matrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<Base<F>> 
SafeHPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    DistMatrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return det::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return det::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<Base<F>> 
SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<Base<F>> 
SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
inline F Determinant( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline Base<F> HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline Base<F> 
HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline Base<F> 
HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

} // namespace elem

#endif // ifndef ELEM_DETERMINANT_HPP
