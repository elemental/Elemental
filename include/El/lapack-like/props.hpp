/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PROPS_HPP
#define EL_PROPS_HPP

namespace El {

// Condition number
// ================
template<typename F>
Base<F> Condition( const Matrix<F>& A, NormType type=TWO_NORM );
template<typename F,Dist U,Dist V>
Base<F> Condition( const DistMatrix<F,U,V>& A, NormType type=TWO_NORM );

// Determinant
// ===========
template<typename F>
SafeProduct<F> SafeDeterminant( const Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( const DistMatrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<F> SafeDeterminant( DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminanti
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const DistMatrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
F Determinant( const Matrix<F>& A );
template<typename F>
F Determinant( const DistMatrix<F>& A );
template<typename F>
F Determinant( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
F Determinant( DistMatrix<F>& A, bool canOverwrite=false );

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const DistMatrix<F>& A );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false );

// Inertia
// =======
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT );
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType=BUNCH_PARLETT );

// Norm
// ====
template<typename F>
Base<F> Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> Norm( const DistMatrix<F,U,V>& A, NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> SymmetricNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> SymmetricNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, 
  NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> HermitianNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F,Dist U,Dist V>
Base<F> HermitianNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, 
  NormType type=FROBENIUS_NORM );

// Entrywise norm
// --------------
template<typename F>
Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p );

// Entrywise one-norm
// ------------------
template<typename F>
Base<F> EntrywiseOneNorm( const Matrix<F>& A );
template<typename F>
Base<F> EntrywiseOneNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianEntrywiseOneNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianEntrywiseOneNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricEntrywiseOneNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricEntrywiseOneNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Frobenius norm
// --------------
template<typename F>
Base<F> FrobeniusNorm( const Matrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Infinity norm
// -------------
template<typename F>
Base<F> InfinityNorm( const Matrix<F>& A );
template<typename F>
Base<F> InfinityNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianInfinityNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricInfinityNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Ky-Fan norms
// ------------
template<typename F>
Base<F> KyFanNorm( const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> KyFanNorm( const DistMatrix<F,U,V>& A, Int k );

template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k );

template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F,Dist U,Dist V>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k );

// Max norm
// --------
template<typename F>
Base<F> MaxNorm( const Matrix<F>& A );
template<typename F>
Base<F> MaxNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianMaxNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Nuclear norm
// ------------
template<typename F>
Base<F> NuclearNorm( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> NuclearNorm( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

// One norm
// --------
template<typename F>
Base<F> OneNorm( const Matrix<F>& A );
template<typename F>
Base<F> OneNorm( const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> HermitianOneNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianOneNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> SymmetricOneNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricOneNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

// Schatten norm
// -------------
template<typename F>
Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> SchattenNorm( const DistMatrix<F,U,V>& A, Base<F> p );

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p );

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F,Dist U,Dist V>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p );

// Two-norm estimate
// -----------------
template<typename F>
Base<F> TwoNormEstimate
( const Matrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> TwoNormEstimate
( const DistMatrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );

template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const DistMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );

template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const DistMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );

// Two norm
// --------
template<typename F>
Base<F> TwoNorm( const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> TwoNorm( const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> HermitianTworNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> HermitianTworNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F,Dist U,Dist V>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

// Zero "norm"
// -----------
template<typename F>
Int ZeroNorm( const Matrix<F>& A );
template<typename F>
Int ZeroNorm( const AbstractDistMatrix<F>& A );

// Trace
// =====
template<typename F>
F Trace( const Matrix<F>& A );
template<typename F>
F Trace( const DistMatrix<F>& A );

} // namespace El

#include "./props/impl.hpp"

#endif // ifndef EL_PROPS_HPP
