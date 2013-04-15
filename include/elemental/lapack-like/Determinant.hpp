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
    PushCallStack("SafeDeterminant");
#endif
    Matrix<F> B( A );
    SafeProduct<F> det = determinant::LUPartialPiv( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SafeDeterminant");
#endif
    DistMatrix<F> B( A );
    SafeProduct<F> det = determinant::LUPartialPiv( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    Matrix<F> B( A );
    SafeProduct<F> det = hpd_determinant::Cholesky( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A );
    SafeProduct<F> det = hpd_determinant::Cholesky( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("SafeDeterminant");
#endif
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    SafeProduct<F> det = determinant::LUPartialPiv( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("SafeDeterminant");
#endif
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    SafeProduct<F> det = determinant::LUPartialPiv( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    Matrix<F> B;
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    SafeProduct<F> det = hpd_determinant::Cholesky( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline SafeProduct<F> 
SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("SafeHPDDeterminant");
#endif
    DistMatrix<F> B( A.Grid() );
    if( canOverwrite )
        View( B, A );
    else
        B = A;
    SafeProduct<F> det = hpd_determinant::Cholesky( uplo, B ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline F Determinant( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Determinant");
#endif
    SafeProduct<F> safeDet = Determinant( A );
    F det = safeDet.rho * Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline F Determinant( const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Determinant");
#endif
    SafeProduct<F> safeDet = Determinant( A );
    F det = safeDet.rho * Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type HPDDeterminant
( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = HPDDeterminant( uplo, A );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type HPDDeterminant
( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = HPDDeterminant( uplo, A );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline F Determinant( Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("Determinant");
#endif
    SafeProduct<F> safeDet = Determinant( A, canOverwrite );
    F det = safeDet.rho * Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline F Determinant( DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("Determinant");
#endif
    SafeProduct<F> safeDet = Determinant( A, canOverwrite );
    F det = safeDet.rho * Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type 
HPDDeterminant( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = HPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

template<typename F>
inline typename Base<F>::type 
HPDDeterminant( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
#ifndef RELEASE
    PushCallStack("HPDDeterminant");
#endif
    SafeProduct<F> safeDet = HPDDeterminant( uplo, A, canOverwrite );
    typename Base<F>::type det = Exp(safeDet.kappa*safeDet.n);
#ifndef RELEASE
    PopCallStack();
#endif
    return det;
}

} // namespace elem

#endif // ifndef LAPACK_DETERMINANT_HPP
