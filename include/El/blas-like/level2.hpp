/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS2_HPP
#define EL_BLAS2_HPP

namespace El {

// Gemv
// ====
template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, Matrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
                 DistMatrix<T>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
  T beta,        DistMatrix<T,VC,STAR>& y );

template<typename T>
void Gemv
( Orientation orientation,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T,VC,STAR>& x,
                 DistMatrix<T,VC,STAR>& y );

template<typename T,Dist AColDist,Dist ARowDist,
                    Dist xColDist,Dist xRowDist,
                    Dist yColDist,Dist yRowDist>
inline void LocalGemv
( Orientation orientation,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A,
           const DistMatrix<T,xColDist,xRowDist>& x,
  T beta,        DistMatrix<T,yColDist,yRowDist>& y )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGemv"))
    // TODO: Add error checking here
    Gemv
    ( orientation ,
      alpha, A.LockedMatrix(), x.LockedMatrix(),
      beta,                    y.Matrix() );
}

// Ger
// ===
template<typename T>
void Ger( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

template<typename T>
void Ger
( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A );

template<typename T,Dist xColDist,Dist xRowDist,
                    Dist yColDist,Dist yRowDist,
                    Dist AColDist,Dist ARowDist>
inline void LocalGer
( T alpha, const DistMatrix<T,xColDist,xRowDist>& x,
           const DistMatrix<T,yColDist,yRowDist>& y,
                 DistMatrix<T,AColDist,ARowDist>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalGer"))
    // TODO: Add error checking here
    Ger( alpha, x.LockedMatrix(), y.LockedMatrix(), A.Matrix() );
}

// Geru
// ====
template<typename T>
void Geru( T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

template<typename T>
void Geru
( T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A );

// Hemv
// ====
template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y );

template<typename T>
void Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y );

// Her2
// ====
template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A );

template<typename T>
void Her2
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A );

// Her
// ===
template<typename T>
void Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A );

template<typename T>
void Her
( UpperOrLower uplo, T alpha, const DistMatrix<T>& x, DistMatrix<T>& A );

// QuasiTrsv
// =========
template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& x,
  bool checkIfSingular=false );

template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation,
  const DistMatrix<F>& A, DistMatrix<F>& x, bool checkIfSingular=false );

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
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y, bool conjugate=false );

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
  bool conjugate=false );

template<typename T>
void LocalRowAccumulate
( UpperOrLower uplo, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false );

} // namespace symv

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
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y,
                 DistMatrix<T>& A, bool conjugate=false );

// Syr
// ===
template<typename T>
void Syr
( UpperOrLower uplo, 
  T alpha, const Matrix<T>& x, Matrix<T>& A, bool conjugate=false );

template<typename T>
void Syr
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x, DistMatrix<T>& A, bool conjugate=false );

// Trmv
// ====
template<typename T>
void Trmv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<T>& A, Matrix<T>& x );
// TODO: Implement distributed version

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
  T alpha, const DistMatrix<T>& X, const DistMatrix<T>& Y, DistMatrix<T>& A,
  bool conjugate=false );

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
  T alpha, const DistMatrix<T>& x, const DistMatrix<T>& y, DistMatrix<T>& A,
  bool conjugate=false );

// Trsv
// ====
template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const Matrix<F>& A, Matrix<F>& x );

template<typename F>
void Trsv
( UpperOrLower uplo, Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F>& A, DistMatrix<F>& x );

} // namespace El

#endif // ifndef EL_BLAS2_HPP
