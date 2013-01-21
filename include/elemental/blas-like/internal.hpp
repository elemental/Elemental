/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_INTERNAL_HPP
#define BLAS_INTERNAL_HPP

namespace elem {
namespace internal {

//----------------------------------------------------------------------------//
// Local BLAS-like: Level 2                                                   //
//----------------------------------------------------------------------------//

template<typename T,Distribution AColDist,Distribution ARowDist,
                    Distribution xColDist,Distribution xRowDist,
                    Distribution yColDist,Distribution yRowDist>
void LocalGemv
( Orientation orientation, 
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A, 
           const DistMatrix<T,xColDist,xRowDist>& x,
  T beta,        DistMatrix<T,yColDist,yRowDist>& y );

template<typename T,Distribution xColDist,Distribution xRowDist,
                    Distribution yColDist,Distribution yRowDist,
                    Distribution AColDist,Distribution ARowDist>
inline void LocalGer
( T alpha, const DistMatrix<T,xColDist,xRowDist>& x, 
           const DistMatrix<T,yColDist,yRowDist>& y,
                 DistMatrix<T,AColDist,ARowDist>& A );


//----------------------------------------------------------------------------//
// Local BLAS-like: Level 3                                                   //
//----------------------------------------------------------------------------//

template<typename T,Distribution AColDist,Distribution ARowDist,
                    Distribution BColDist,Distribution BRowDist,
                    Distribution CColDist,Distribution CRowDist>
void LocalGemm
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,AColDist,ARowDist>& A, 
           const DistMatrix<T,BColDist,BRowDist>& B,
  T beta,        DistMatrix<T,CColDist,CRowDist>& C );

template<typename T>
void LocalTrtrmm
( Orientation orientation, UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A );

template<typename T>
void LocalTrdtrmm
( Orientation orientation, UpperOrLower uplo, DistMatrix<T,STAR,STAR>& A );

template<typename T,Distribution BColDist,Distribution BRowDist>
void LocalTrmm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  T alpha, const DistMatrix<T,STAR,STAR>& A,
                 DistMatrix<T,BColDist,BRowDist>& B );

template<typename F,Distribution XColDist,Distribution XRowDist>
void LocalTrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A, 
                 DistMatrix<F,XColDist,XRowDist>& X,
  bool checkIfSingular=false );

template<typename F>
void LocalTrtrsm
( LeftOrRight side, UpperOrLower uplo, 
  Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,STAR>& A, 
                 DistMatrix<F,STAR,STAR>& X,
  bool checkIfSingular=false );

// TODO: Finish adding wrappers for Local BLAS-like routines

//----------------------------------------------------------------------------//
// Distributed BLAS-like helpers: Level 2                                     //
//----------------------------------------------------------------------------//
            
// This is for the case where x is a column vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T x[MC,* ]
template<typename T>
void LocalSymvColAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate=false );

// This is for the case where x is a column vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trius(A)[MC,MR])^T x[MC,* ]
template<typename T>
void LocalSymvColAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate=false );

// This is for the case where x is a row vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T (x[* ,MC])^T
template<typename T>
void LocalSymvRowAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false );

// This is for the case where x is a row vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trius(A)[MC,MR])^T (x[* ,MC])^T
template<typename T>
void LocalSymvRowAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false );

//----------------------------------------------------------------------------//
// Distributed BLAS-like helpers: Level 3                                     //
//----------------------------------------------------------------------------//

template<typename T>
void LocalSymmetricAccumulateLL
( Orientation orientation, T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BHermOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR
);

template<typename T>
void LocalSymmetricAccumulateLU
( Orientation orientation, T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,  STAR>& B_MC_STAR,
  const DistMatrix<T,STAR,MR  >& BHermOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR
);

template<typename T>
void LocalSymmetricAccumulateRL
( Orientation orientation, T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BHermOrTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZHermOrTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZHermOrTrans_MR_STAR
);

template<typename T>
void LocalSymmetricAccumulateRU
( Orientation orientation, T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BHermOrTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZHermOrTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZHermOrTrans_MR_STAR
);

template<typename T>
void LocalTrmmAccumulateLLN
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& L,
  const DistMatrix<T,STAR,MR  >& XHermOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR );

template<typename T>
void LocalTrmmAccumulateLLT
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& L,
  const DistMatrix<T,MC,  STAR>& X_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR );

template<typename T>
void LocalTrmmAccumulateLUN
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& U,
  const DistMatrix<T,STAR,MR  >& XHermOrTrans_STAR_MR,
        DistMatrix<T,MC,  STAR>& Z_MC_STAR );

template<typename T>
void LocalTrmmAccumulateLUT
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& U,
  const DistMatrix<T,MC,  STAR>& X_MC_STAR,
        DistMatrix<T,MR,  STAR>& Z_MR_STAR );

template<typename T>
void LocalTrmmAccumulateRLN
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& L,
  const DistMatrix<T,STAR,MC  >& X_STAR_MC,
        DistMatrix<T,MR,  STAR>& ZHermOrTrans_MR_STAR );

template<typename T>
void LocalTrmmAccumulateRLT
( UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& L,
  const DistMatrix<T,MR,  STAR>& XHermOrTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZHermOrTrans_MC_STAR );

template<typename T>
void LocalTrmmAccumulateRUN
( Orientation orientation, UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& U,
  const DistMatrix<T,STAR,MC  >& X_STAR_MC,
        DistMatrix<T,MR,  STAR>& ZHermOrTrans_MR_STAR );

template<typename T>
void LocalTrmmAccumulateRUT
( UnitOrNonUnit diag, T alpha, 
  const DistMatrix<T>& U,
  const DistMatrix<T,MR,  STAR>& XHermOrTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZHermOrTrans_MC_STAR );

// Triangular rank-k Update:
// tril(C) := alpha tril( A B ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A B ) + beta triu(C)

template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T>& C );

// Triangular rank-k Update:
// tril(C) := alpha tril( A B^{T/H} ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A B^{T/H} ) + beta triu(C)

template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C );

// Triangular rank-k Update:
// tril(C) := alpha tril( A^{T/H} B ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A^{T/H} B ) + beta triu(C)

template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
  T beta,        DistMatrix<T>& C );

// Triangular rank-k Update:
// tril(C) := alpha tril( A^{T/H} B^{T/H} ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A^{T/H} B^{T/H} ) + beta triu(C)

template<typename T>
void LocalTrrk
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
  T beta,        DistMatrix<T>& C );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B + C D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B + C D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B, 
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B + C D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B + C D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B + C^{T/H} D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B + C^{T/H} D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B + C^{T/H} D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B + C^{T/H} D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B^{T/H} + C D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B^{T/H} + C D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B^{T/H} + C D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B^{T/H} + C D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B^{T/H} + C^{T/H} D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B^{T/H} + C^{T/H} D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A B^{T/H} + C^{T/H} D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A B^{T/H} + C^{T/H} D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,STAR>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B + C D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B + C D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B + C D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B + C D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC  >& A, const DistMatrix<T,STAR,MR  >& B,
           const DistMatrix<T,MC,  STAR>& C, const DistMatrix<T,MR,  STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B + C^{T/H} D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B + C^{T/H} D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B + C^{T/H} D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B + C^{T/H} D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,STAR,MR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B^{T/H} + C D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B^{T/H} + C D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B^{T/H} + C D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B^{T/H} + C D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,MC,STAR>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B^{T/H} + C^{T/H} D ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B^{T/H} + C^{T/H} D ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B,
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,STAR,MR>& D,
  T beta,        DistMatrix<T>& E );

// Triangular rank-2k Update:
// tril(E) := alpha tril( A^{T/H} B^{T/H} + C^{T/H} D^{T/H} ) + beta tril(E)
//   or
// triu(E) := alpha triu( A^{T/H} B^{T/H} + C^{T/H} D^{T/H} ) + beta triu(E)

template<typename T>
void LocalTrr2k
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfB,
  Orientation orientationOfC,
  Orientation orientationOfD,
  T alpha, const DistMatrix<T,STAR,MC>& A, const DistMatrix<T,MR,STAR>& B, 
           const DistMatrix<T,STAR,MC>& C, const DistMatrix<T,MR,STAR>& D,
  T beta,        DistMatrix<T>& E );

// Left, Lower, Normal Trsm
template<typename F>
void TrsmLLNLarge
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLLNMedium
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLLNSmall
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular=false );

// Left, Lower, (Conjugate)Transpose Trsm
template<typename F>
void TrsmLLTLarge
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLLTMedium
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& L, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLLTSmall
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,VC,STAR>& L, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLLTSmall
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,VR>& L, DistMatrix<F,VR,STAR>& X,
  bool checkIfSingular=false );

// Left, Upper, Normal Trsm
template<typename F>
void TrsmLUNLarge
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLUNMedium
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLUNSmall
( UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,VC,STAR>& U, DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular=false );

// Left, Upper, (Conjugate)Transpose Trsm
template<typename F>
void TrsmLUTLarge
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLUTMedium
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular=false );
template<typename F>
void TrsmLUTSmall
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,VR>& U, DistMatrix<F,VR,STAR>& X,
  bool checkIfSingular=false );

} // internal
} // elem

#endif // ifndef BLAS_INTERNAL_HPP
