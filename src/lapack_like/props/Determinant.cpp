/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Determinant/Cholesky.hpp"
#include "./Determinant/LUPartialPiv.hpp"

namespace El {

template<typename Field>
SafeProduct<Field> SafeDeterminant( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Field> B( A );
    return det::LUPartialPiv( B );
}

template<typename Field>
SafeProduct<Field> SafeDeterminant( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Field> B( A );
    return det::LUPartialPiv( B );
}

template<typename Field>
SafeProduct<Base<Field>> SafeHPDDeterminant
( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    Matrix<Field> B( A );
    return hpd_det::Cholesky( uplo, B );
}

template<typename Field>
SafeProduct<Base<Field>> SafeHPDDeterminant
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    DistMatrix<Field> B( A );
    return hpd_det::Cholesky( uplo, B );
}

template<typename Field>
SafeProduct<Field> SafeDeterminant( Matrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    Matrix<Field> B;
    if( canOverwrite )
    {
        return det::LUPartialPiv( A );
    }
    else
    {
        Matrix<Field> B( A );
        return det::LUPartialPiv( B );
    }
}

template<typename Field>
SafeProduct<Field>
SafeDeterminant( AbstractDistMatrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    if( canOverwrite )
    {
        return det::LUPartialPiv( A );
    }
    else
    {
        DistMatrix<Field> B( A );
        return det::LUPartialPiv( B );
    }
}

template<typename Field>
SafeProduct<Base<Field>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    if( canOverwrite )
    {
        return hpd_det::Cholesky( uplo, A );
    }
    else
    {
        Matrix<Field> B( A );
        return hpd_det::Cholesky( uplo, B );
    }
}

template<typename Field>
SafeProduct<Base<Field>> SafeHPDDeterminant
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    if( canOverwrite )
    {
        return hpd_det::Cholesky( uplo, A );
    }
    else
    {
        DistMatrix<Field> B( A );
        return hpd_det::Cholesky( uplo, B );
    }
}

template<typename Field>
Field Determinant( const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Field> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Field Determinant( const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Field> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Base<Field> HPDDeterminant( UpperOrLower uplo, const Matrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Base<Field>
HPDDeterminant( UpperOrLower uplo, const AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Field Determinant( Matrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Field> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Field Determinant( AbstractDistMatrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Field> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Base<Field> HPDDeterminant
( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet =
      SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename Field>
Base<Field> HPDDeterminant
( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite )
{
    EL_DEBUG_CSE
    SafeProduct<Base<Field>> safeDet =
      SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

#define PROTO(Field) \
  template SafeProduct<Field> SafeDeterminant( const Matrix<Field>& A ); \
  template SafeProduct<Field> SafeDeterminant \
  ( const AbstractDistMatrix<Field>& A ); \
  template SafeProduct<Field> SafeDeterminant \
  ( Matrix<Field>& A, bool canOverwrite ); \
  template SafeProduct<Field> SafeDeterminant \
  ( AbstractDistMatrix<Field>& A, bool canOverwrite ); \
  template SafeProduct<Base<Field>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template SafeProduct<Base<Field>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template SafeProduct<Base<Field>> SafeHPDDeterminant \
  ( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite ); \
  template SafeProduct<Base<Field>> SafeHPDDeterminant \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite ); \
  template Field Determinant( const Matrix<Field>& A ); \
  template Field Determinant( const AbstractDistMatrix<Field>& A ); \
  template Field Determinant( Matrix<Field>& A, bool canOverwrite ); \
  template Field Determinant \
  ( AbstractDistMatrix<Field>& A, bool canOverwrite ); \
  template Base<Field> HPDDeterminant \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template Base<Field> HPDDeterminant \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template Base<Field> HPDDeterminant \
  ( UpperOrLower uplo, Matrix<Field>& A, bool canOverwrite ); \
  template Base<Field> HPDDeterminant \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, bool canOverwrite ); \
  template SafeProduct<Base<Field>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const Matrix<Field>& A ); \
  template SafeProduct<Base<Field>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A ); \
  template SafeProduct<Field> det::AfterLUPartialPiv \
  ( const Matrix<Field>& A, const Permutation& P ); \
  template SafeProduct<Field> det::AfterLUPartialPiv \
  ( const AbstractDistMatrix<Field>& A, const DistPermutation& P );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
