/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_DETERMINANT_HPP
#define LAPACK_DETERMINANT_HPP

#include "elemental/lapack-like/Determinant/Cholesky.hpp"
#include "elemental/lapack-like/Determinant/LUPartialPiv.hpp"

namespace elem {

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SafeDeterminant");
#endif
    Matrix<F> B( A );
    return determinant::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SafeDeterminant");
#endif
    DistMatrix<F> B( A );
    return determinant::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SafeHPDDeterminant");
#endif
    Matrix<F> B( A );
    return hpd_determinant::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A );
    return hpd_determinant::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("SafeDeterminant");
#endif
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return determinant::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("SafeDeterminant");
#endif
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return determinant::LUPartialPiv( B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("SafeHPDDeterminant");
#endif
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return hpd_determinant::Cholesky( uplo, B ); 
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    return hpd_determinant::Cholesky( uplo, B ); 
}

template<typename F>
inline F Determinant( const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Determinant");
#endif
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Determinant");
#endif
    SafeProduct<F> safeDet = SafeDeterminant( A );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline BASE(F) HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline BASE(F) HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("Determinant");
#endif
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline F Determinant( DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("Determinant");
#endif
    SafeProduct<F> safeDet = SafeDeterminant( A, canOverwrite );
    return safeDet.rho * Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline BASE(F) 
HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

template<typename F>
inline BASE(F) 
HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    CallStackEntry entry("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return Exp(safeDet.kappa*safeDet.n);
}

} // namespace elem

#endif // ifndef LAPACK_DETERMINANT_HPP
