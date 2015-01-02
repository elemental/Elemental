/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Determinant/Cholesky.hpp"
#include "./Determinant/LUPartialPiv.hpp"

namespace El {

template<typename F>
SafeProduct<F>  SafeDeterminant( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    Matrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

// TODO: Switch to AbstractDistMatrix
template<typename F>
SafeProduct<F> SafeDeterminant( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    DistMatrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    Matrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    DistMatrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
SafeProduct<F> SafeDeterminant( Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    Matrix<F> B;
    if( canOverwrite )
    {
        return det::LUPartialPiv( A ); 
    }
    else
    {
        Matrix<F> B( A );
        return det::LUPartialPiv( B ); 
    }
}

template<typename F>
SafeProduct<F> SafeDeterminant( AbstractDistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("SafeDeterminant"))
    if( canOverwrite )
    {
        return det::LUPartialPiv( A ); 
    }
    else
    {
        DistMatrix<F> B( A );
        return det::LUPartialPiv( B ); 
    }
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    if( canOverwrite )
    {
        return hpd_det::Cholesky( uplo, A ); 
    }
    else
    {
        Matrix<F> B( A );
        return hpd_det::Cholesky( uplo, B ); 
    }
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("SafeHPDDeterminant"))
    if( canOverwrite )
    {
        return hpd_det::Cholesky( uplo, A ); 
    }
    else
    {
        DistMatrix<F> B( A );
        return hpd_det::Cholesky( uplo, B ); 
    }
}

template<typename F>
F Determinant( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( AbstractDistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

#define PROTO(F) \
  template SafeProduct<F> SafeDeterminant( const Matrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant( const AbstractDistMatrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant \
  ( Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<F> SafeDeterminant \
  ( AbstractDistMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite ); \
  \
  template F Determinant( const Matrix<F>& A ); \
  template F Determinant( const AbstractDistMatrix<F>& A ); \
  template F Determinant( Matrix<F>& A, bool canOverwrite ); \
  template F Determinant( AbstractDistMatrix<F>& A, bool canOverwrite ); \
  \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const Matrix<F>& A, const Matrix<Int>& p ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<Int>& p );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
