/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PERM_HPP
#define EL_PERM_HPP

namespace El {

struct PermutationMeta
{
    Int align;
    mpi::Comm comm;

    // Will treat vector lengths as one
    vector<int> sendCounts, sendDispls,
                recvCounts, recvDispls;

    vector<int> sendIdx, sendRanks,
                recvIdx, recvRanks;

    int TotalSend() const { return sendCounts.back()+sendDispls.back(); }
    int TotalRecv() const { return recvCounts.back()+recvDispls.back(); }

    void ScaleUp( Int length )
    {
        const int p = sendCounts.size();
        for( int q=0; q<p; ++q )
        {
            sendCounts[q] *= length;
            sendDispls[q] *= length;
            recvCounts[q] *= length;
            recvDispls[q] *= length;
        }
    }
    void ScaleDown( Int length )
    {
        const int p = sendCounts.size();
        for( int q=0; q<p; ++q )
        {
            sendCounts[q] /= length;
            sendDispls[q] /= length;
            recvCounts[q] /= length;
            recvDispls[q] /= length;
        }
    }

    PermutationMeta()
    : align(0), comm(mpi::COMM_SELF), 
      sendCounts(1,0), sendDispls(1,0), 
      recvCounts(1,0), recvDispls(1,0)
    { }

    PermutationMeta
    ( const AbstractDistMatrix<Int>& p,
      const AbstractDistMatrix<Int>& pInv );
};

// Apply column pivots
// ===================
template<typename T>
void ApplyColPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T>
void ApplyColPivots
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& pivots, 
  Int offset=0 );

template<typename T>
void ApplyInverseColPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T>
void ApplyInverseColPivots
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& pivots, 
  Int offset=0 );

// Apply row pivots
// ================
template<typename T>
void ApplyRowPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T>
void ApplyRowPivots
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& pivots,
  Int offset=0 );

template<typename T>
void ApplyInverseRowPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T>
void ApplyInverseRowPivots
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& pivots, 
  Int offset=0 );

// Apply symmetric pivots
// ======================
template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate=false, 
  Int offset=0 );
template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo, AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots, bool conjugate=false, 
  Int offset=0 );

template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate=false, 
  Int offset=0 );
template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots, bool conjugate=false, 
  Int offset=0 );

// Explicit permutation
// ====================
void ExplicitPermutation
( const Matrix<Int>& p, Matrix<Int>& P );
void ExplicitPermutation
( const AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& P );

// Invert permutation
// ==================
void InvertPermutation( const Matrix<Int>& p, Matrix<Int>& pInv );
void InvertPermutation
( const AbstractDistMatrix<Int>& p, AbstractDistMatrix<Int>& pInv );

// Parity of a permutation
// =======================
bool PermutationParity( const Matrix<Int>& p );
bool PermutationParity( const AbstractDistMatrix<Int>& p );

// Permute columns
// ===============
template<typename T>
void PermuteCols( Matrix<T>& A, const Matrix<Int>& p );
template<typename T>
void PermuteCols
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& p );

template<typename T>
void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& p );
template<typename T>
void InversePermuteCols
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& p );

template<typename T>
void PermuteCols
( Matrix<T>& A, const Matrix<Int>& p, const Matrix<Int>& pInv );
template<typename T>
void PermuteCols
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p,
  const AbstractDistMatrix<Int>& pInv );

template<typename T>
void PermuteCols( AbstractDistMatrix<T>& A, const PermutationMeta& meta );

// Permute rows
// ============
template<typename T>
void PermuteRows( Matrix<T>& A, const Matrix<Int>& p );
template<typename T>
void PermuteRows
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& p );

template<typename T>
void InversePermuteRows( Matrix<T>& A, const Matrix<Int>& p );
template<typename T>
void InversePermuteRows
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& p );

template<typename T>
void PermuteRows
( Matrix<T>& A, const Matrix<Int>& p, const Matrix<Int>& pInv );
template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p,
  const AbstractDistMatrix<Int>& pInv );

template<typename T>
void PermuteRows( AbstractDistMatrix<T>& A, const PermutationMeta& meta );

// Parity of a sequence of pivots
// ==============================
bool PivotParity( const Matrix<Int>& p, Int pivotOffset=0 );
bool PivotParity( const AbstractDistMatrix<Int>& p, Int pivotOffset=0 );

// Convert a pivot sequence to a partial permutation vector
// ========================================================
void PivotsToPartialPermutation
( const Matrix<Int>& pivots, Matrix<Int>& p, Matrix<Int>& pInv,
  Int offset=0 );
void PivotsToPartialPermutation
( const AbstractDistMatrix<Int>& pivots,
        AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv, Int offset=0 );

// Convert a pivot sequence to a permutation vector
// ================================================
void PivotsToPermutation
( const Matrix<Int>& pivots, Matrix<Int>& p, Int offset=0 );
void PivotsToPermutation
( const AbstractDistMatrix<Int>& pivots, AbstractDistMatrix<Int>& p,
  Int offset=0 );

void PivotsToInversePermutation
( const Matrix<Int>& pivots, Matrix<Int>& pInv, Int offset=0 );
void PivotsToInversePermutation
( const AbstractDistMatrix<Int>& pivots, AbstractDistMatrix<Int>& pInv,
  Int offset=0 );

} // namespace El

#endif // ifndef EL_PERM_HPP
