/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS2_HPP
#define EL_BLAS2_HPP

namespace El {

// Gemv
// ====
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& x,
  T beta,        Matrix<T>& y );
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A,
           const Matrix<T>& x,
                 Matrix<T>& y );
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& x,
  T beta,        ElementalMatrix<T>& y );
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const ElementalMatrix<T>& A,
           const ElementalMatrix<T>& x,
                 ElementalMatrix<T>& y );
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR,BLOCK>& A,
           const DistMatrix<T,MC,MR,BLOCK>& x,
  T beta,        DistMatrix<T,MC,MR,BLOCK>& y );
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR,BLOCK>& A,
           const DistMatrix<T,MC,MR,BLOCK>& x,
                 DistMatrix<T,MC,MR,BLOCK>& y );
template<typename T>
void LocalGemv
( Orientation orientation,
  T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& x,
  T beta,        ElementalMatrix<T>& y );

// Ger
// ===
template<typename T>
void Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );
template<typename T>
void Ger
( T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y, 
                 ElementalMatrix<T>& A );
template<typename T>
void LocalGer
( T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y,
                 ElementalMatrix<T>& A );

// Geru
// ====
template<typename T>
void Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

template<typename T>
void Geru
( T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y, 
                 ElementalMatrix<T>& A );

// Hemv
// ====
template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& x,
  T beta,        ElementalMatrix<T>& y,
  const SymvCtrl<T>& ctrl=SymvCtrl<T>() );

// Her
// ===
template<typename T>
void Her( UpperOrLower uplo, Base<T> alpha, const Matrix<T>& x, Matrix<T>& A );

template<typename T>
void Her
( UpperOrLower uplo, 
  Base<T> alpha, const ElementalMatrix<T>& x, ElementalMatrix<T>& A );

// Her2
// ====
template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y,
                 ElementalMatrix<T>& A );

// QuasiTrsv
// =========
template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& x,
  bool checkIfSingular=false );

template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A, ElementalMatrix<F>& x, 
  bool checkIfSingular=false );

// Symv
// ====
template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y,
  bool conjugate=false );

template<typename T>
void Symv
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& x,
  T beta,        ElementalMatrix<T>& y, bool conjugate=false, 
  const SymvCtrl<T>& ctrl=SymvCtrl<T>() );

// namespace symv
// --------------
namespace symv {

template<typename T>
void LocalColAccumulate
( UpperOrLower uplo, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate=false, const SymvCtrl<T>& ctrl=SymvCtrl<T>() );

template<typename T>
void LocalRowAccumulate
( UpperOrLower uplo, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false, const SymvCtrl<T>& ctrl=SymvCtrl<T>() );

} // namespace symv

// Syr
// ===
template<typename T>
void Syr
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, Matrix<T>& A, bool conjugate=false );

template<typename T>
void Syr
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& x, ElementalMatrix<T>& A, 
  bool conjugate=false );

// Syr2
// ====
template<typename T>
void Syr2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A,
  bool conjugate=false );

template<typename T>
void Syr2
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y,
                 ElementalMatrix<T>& A, bool conjugate=false );

// Trmv
// ====
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<T>& A, Matrix<T>& x );
// TODO: Implement distributed version

// Trr
// ===
// A := A + alpha x y'
template<typename T>
void Trr
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A,
  bool conjugate=false );

template<typename T>
void Trr
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& x, const ElementalMatrix<T>& y, 
  ElementalMatrix<T>& A, bool conjugate=false );

// Trr2
// ====
// A := A + alpha X Y'
template<typename T>
void Trr2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& X, const Matrix<T>& Y, Matrix<T>& A,
  bool conjugate=false );

template<typename T>
void Trr2
( UpperOrLower uplo,
  T alpha, const ElementalMatrix<T>& X, const ElementalMatrix<T>& Y, 
  ElementalMatrix<T>& A, bool conjugate=false );

// Trsv
// ====
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x );

template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& x );

} // namespace El

#endif // ifndef EL_BLAS2_HPP
