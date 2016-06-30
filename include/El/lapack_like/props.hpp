/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PROPS_HPP
#define EL_PROPS_HPP

namespace El {

// Condition number
// ================
template<typename F>
Base<F> Condition( const Matrix<F>& A, NormType type=TWO_NORM );
template<typename F>
Base<F> Condition( const ElementalMatrix<F>& A, NormType type=TWO_NORM );

template<typename F>
Base<F> FrobeniusCondition( const Matrix<F>& A );
template<typename F>
Base<F> FrobeniusCondition( const ElementalMatrix<F>& A );

template<typename F>
Base<F> InfinityCondition( const Matrix<F>& A );
template<typename F>
Base<F> InfinityCondition( const ElementalMatrix<F>& A );

template<typename F>
Base<F> MaxCondition( const Matrix<F>& A );
template<typename F>
Base<F> MaxCondition( const ElementalMatrix<F>& A );

template<typename F>
Base<F> OneCondition( const Matrix<F>& A );
template<typename F>
Base<F> OneCondition( const ElementalMatrix<F>& A );

template<typename F>
Base<F> TwoCondition( const Matrix<F>& A );
template<typename F>
Base<F> TwoCondition( const ElementalMatrix<F>& A );

// Determinant
// ===========
template<typename F>
SafeProduct<F> SafeDeterminant( const Matrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant( const ElementalMatrix<F>& A );
template<typename F>
SafeProduct<F> SafeDeterminant
( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<F> SafeDeterminant
( ElementalMatrix<F>& A, bool canOverwrite=false );

template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, const ElementalMatrix<F>& A );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
SafeProduct<Base<F>> SafeHPDDeterminant
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite=false );

template<typename F>
F Determinant( const Matrix<F>& A );
template<typename F>
F Determinant( const ElementalMatrix<F>& A );
template<typename F>
F Determinant( Matrix<F>& A, bool canOverwrite=false );
template<typename F>
F Determinant( ElementalMatrix<F>& A, bool canOverwrite=false );

template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HPDDeterminant( UpperOrLower uplo, const ElementalMatrix<F>& A );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
Base<F> HPDDeterminant
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite=false );

namespace hpd_det {

template<typename F>
SafeProduct<Base<F>> AfterCholesky
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
SafeProduct<Base<F>> AfterCholesky
( UpperOrLower uplo, const ElementalMatrix<F>& A );

} // namespace hpd_det

namespace det {

template<typename F>
SafeProduct<F> AfterLUPartialPiv
( const Matrix<F>& A, const Permutation& P );
template<typename F>
SafeProduct<F> AfterLUPartialPiv
( const ElementalMatrix<F>& A, const DistPermutation& P );

} // namespace det

// Inertia
// =======
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, Matrix<F>& A,
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
InertiaType Inertia
( UpperOrLower uplo, ElementalMatrix<F>& A, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

// Norm
// ====
template<typename F>
Base<F> Norm( const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F>
Base<F> Norm( const ElementalMatrix<F>& A, NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> SymmetricNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F>
Base<F> SymmetricNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, 
  NormType type=FROBENIUS_NORM );

template<typename F>
Base<F> HermitianNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM );
template<typename F>
Base<F> HermitianNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, 
  NormType type=FROBENIUS_NORM );

// Entrywise norm
// --------------
template<typename F>
Base<F> EntrywiseNorm( const Matrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> EntrywiseNorm( const AbstractDistMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> EntrywiseNorm( const SparseMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> EntrywiseNorm( const DistSparseMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> EntrywiseNorm( const DistMultiVec<F>& A, Base<F> p=1 );

template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> HermitianEntrywiseNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p=1 );

template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const SparseMatrix<F>& A, Base<F> p=1 );
template<typename F>
Base<F> SymmetricEntrywiseNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A, Base<F> p=1 );

// Frobenius norm
// --------------
template<typename F>
Base<F> FrobeniusNorm( const Matrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const AbstractDistMatrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const SparseMatrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const DistSparseMatrix<F>& A );
template<typename F>
Base<F> FrobeniusNorm( const DistMultiVec<F>& A );

template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );
template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const SparseMatrix<F>& A );
template<typename F>
Base<F> HermitianFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A );

template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const AbstractDistMatrix<F>& A );
template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const SparseMatrix<F>& A );
template<typename F>
Base<F> SymmetricFrobeniusNorm
( UpperOrLower uplo, const DistSparseMatrix<F>& A );

// Infinity norm
// -------------
template<typename T>
Base<T> InfinityNorm( const Matrix<T>& A );
template<typename T>
Base<T> InfinityNorm( const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> InfinityNorm( const SparseMatrix<T>& A );
template<typename T>
Base<T> InfinityNorm( const DistSparseMatrix<T>& A );

template<typename T>
Base<T> HermitianInfinityNorm
( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> HermitianInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<T>& A );

template<typename T>
Base<T> SymmetricInfinityNorm
( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> SymmetricInfinityNorm
( UpperOrLower uplo, const AbstractDistMatrix<T>& A );

template<typename T>
Base<T> HermitianTridiagInfinityNorm
( const Matrix<Base<T>>& d, const Matrix<T>& e );

// Ky-Fan norms
// ------------
template<typename F>
Base<F> KyFanNorm( const Matrix<F>& A, Int k );
template<typename F>
Base<F> KyFanNorm( const ElementalMatrix<F>& A, Int k );

template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k );

template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k );
template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k );

// Ky-Fan-Schatten norms
// ---------------------
template<typename F>
Base<F> KyFanSchattenNorm( const Matrix<F>& A, Int k, Base<F> p );
template<typename F>
Base<F> KyFanSchattenNorm( const ElementalMatrix<F>& A, Int k, Base<F> p );

template<typename F>
Base<F> HermitianKyFanSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k, Base<F> p );
template<typename F>
Base<F> HermitianKyFanSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k, Base<F> p );

template<typename F>
Base<F> SymmetricKyFanSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Int k, Base<F> p );
template<typename F>
Base<F> SymmetricKyFanSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k, Base<F> p );

// Max norm
// --------
template<typename T>
Base<T> MaxNorm( const Matrix<T>& A );
template<typename T>
Base<T> MaxNorm( const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> MaxNorm( const SparseMatrix<T>& A );
template<typename T>
Base<T> MaxNorm( const DistSparseMatrix<T>& A );
template<typename T>
Base<T> MaxNorm( const DistMultiVec<T>& A );

template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A );
template<typename T>
Base<T> HermitianMaxNorm( UpperOrLower uplo, const DistSparseMatrix<T>& A );

template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const SparseMatrix<T>& A );
template<typename T>
Base<T> SymmetricMaxNorm( UpperOrLower uplo, const DistSparseMatrix<T>& A );

// Nuclear norm
// ------------
template<typename F>
Base<F> NuclearNorm( const Matrix<F>& A );
template<typename F>
Base<F> NuclearNorm( const ElementalMatrix<F>& A );

template<typename F>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianNuclearNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A );

template<typename F>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricNuclearNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A );

// One norm
// --------
template<typename T>
Base<T> OneNorm( const Matrix<T>& A );
template<typename T>
Base<T> OneNorm( const AbstractDistMatrix<T>& A );
template<typename T>
Base<T> OneNorm( const SparseMatrix<T>& A );
template<typename T>
Base<T> OneNorm( const DistSparseMatrix<T>& A );

template<typename T>
Base<T> HermitianOneNorm( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> HermitianOneNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A );

template<typename T>
Base<T> SymmetricOneNorm( UpperOrLower uplo, const Matrix<T>& A );
template<typename T>
Base<T> SymmetricOneNorm( UpperOrLower uplo, const AbstractDistMatrix<T>& A );

template<typename T>
Base<T> HermitianTridiagOneNorm
( const Matrix<Base<T>>& d, const Matrix<T>& e );

// Schatten norm
// -------------
template<typename F>
Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> SchattenNorm( const ElementalMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p );

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p );
template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p );

// Two norm
// --------
template<typename F>
Base<F> TwoNorm( const Matrix<F>& A );
template<typename F>
Base<F> TwoNorm( const ElementalMatrix<F>& A );

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const ElementalMatrix<F>& A );

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const ElementalMatrix<F>& A );

// Zero "norm"
// -----------
template<typename T>
Int ZeroNorm( const Matrix<T>& A, Base<T> tol=0 );
template<typename T>
Int ZeroNorm( const SparseMatrix<T>& A, Base<T> tol=0 );
template<typename T>
Int ZeroNorm( const AbstractDistMatrix<T>& A, Base<T> tol=0 );
template<typename T>
Int ZeroNorm( const DistSparseMatrix<T>& A, Base<T> tol=0 );

// Two-norm estimate
// -----------------
template<typename F>
Base<F> TwoNormEstimate
( const Matrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> TwoNormEstimate
( const ElementalMatrix<F>& A, Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> TwoNormEstimate( const SparseMatrix<F>& A, Int basisSize=15 );
template<typename F>
Base<F> TwoNormEstimate( const DistSparseMatrix<F>& A, Int basisSize=15 );

template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> HermitianTwoNormEstimate
( UpperOrLower uplo, const ElementalMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> HermitianTwoNormEstimate
( const SparseMatrix<F>& A, Int basisSize=15 );
template<typename F>
Base<F> HermitianTwoNormEstimate
( const DistSparseMatrix<F>& A, Int basisSize=15 );

template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const Matrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );
template<typename F>
Base<F> SymmetricTwoNormEstimate
( UpperOrLower uplo, const ElementalMatrix<F>& A, 
  Base<F> tol=1e-6, Int maxIts=1000 );

// Trace
// =====
template<typename T>
T Trace( const Matrix<T>& A );
template<typename T>
T Trace( const AbstractDistMatrix<T>& A );

} // namespace El

#endif // ifndef EL_PROPS_HPP
