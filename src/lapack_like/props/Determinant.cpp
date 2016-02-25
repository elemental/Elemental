/*
   Copyright (c) 2009-2016, Jack Poulson
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
    DEBUG_ONLY(CSE cse("SafeDeterminant"))
    Matrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

template<typename F>
SafeProduct<F> SafeDeterminant( const ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SafeDeterminant"))
    DistMatrix<F> B( A );
    return det::LUPartialPiv( B ); 
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SafeHPDDeterminant"))
    Matrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SafeHPDDeterminant"))
    DistMatrix<F> B( A );
    return hpd_det::Cholesky( uplo, B ); 
}

template<typename F>
SafeProduct<F> SafeDeterminant( Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("SafeDeterminant"))
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
SafeProduct<F> SafeDeterminant( ElementalMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("SafeDeterminant"))
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
    DEBUG_ONLY(CSE cse("SafeHPDDeterminant"))
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
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("SafeHPDDeterminant"))
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
    DEBUG_ONLY(CSE cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( const ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( ElementalMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CSE cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

#define PROTO(F) \
  template SafeProduct<F> SafeDeterminant( const Matrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant( const ElementalMatrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant \
  ( Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<F> SafeDeterminant \
  ( ElementalMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite ); \
  \
  template F Determinant( const Matrix<F>& A ); \
  template F Determinant( const ElementalMatrix<F>& A ); \
  template F Determinant( Matrix<F>& A, bool canOverwrite ); \
  template F Determinant( ElementalMatrix<F>& A, bool canOverwrite ); \
  \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const Matrix<F>& A, const Permutation& P ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const ElementalMatrix<F>& A, const DistPermutation& P );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
