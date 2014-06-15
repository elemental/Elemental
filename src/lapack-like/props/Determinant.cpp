/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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

template<typename F>
SafeProduct<F> SafeDeterminant( const DistMatrix<F>& A )
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
( UpperOrLower uplo, const DistMatrix<F>& A )
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
        View( B, A );
    else
        B = A;
    return det::LUPartialPiv( B ); 
}

template<typename F>
SafeProduct<F> SafeDeterminant( DistMatrix<F>& A, bool canOverwrite )
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
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
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
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite )
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
F Determinant( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Determinant"))
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
F Determinant( const DistMatrix<F>& A )
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
Base<F> HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
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
F Determinant( DistMatrix<F>& A, bool canOverwrite )
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
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("HPDDeterminant"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

#define PROTO(F) \
  template SafeProduct<F> SafeDeterminant( const Matrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant( const DistMatrix<F>& A ); \
  template SafeProduct<F> SafeDeterminant \
  ( Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<F> SafeDeterminant \
  ( DistMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, const DistMatrix<F>& A ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template SafeProduct<Base<F>> SafeHPDDeterminant \
  ( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite ); \
  \
  template F Determinant( const Matrix<F>& A ); \
  template F Determinant( const DistMatrix<F>& A ); \
  template F Determinant( Matrix<F>& A, bool canOverwrite ); \
  template F Determinant( DistMatrix<F>& A, bool canOverwrite ); \
  \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, const DistMatrix<F>& A ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> HPDDeterminant \
  ( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite ); \
  \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template SafeProduct<Base<F>> hpd_det::AfterCholesky \
  ( UpperOrLower uplo, const DistMatrix<F>& A ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const Matrix<F>& A, const Matrix<Int>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,MC,STAR>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,MD,STAR>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,MR,STAR>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,STAR,STAR>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& pPerm ); \
  template SafeProduct<F> det::AfterLUPartialPiv \
  ( const DistMatrix<F>& A, const DistMatrix<Int,VR,STAR>& pPerm );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
