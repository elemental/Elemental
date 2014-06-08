/*
   Copyright (c) 2009-2014, Jack Poulson
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
    std::vector<int> sendCounts, sendDispls,
                     recvCounts, recvDispls;

    std::vector<int> sendIdx, sendRanks,
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

    template<Dist U>
    PermutationMeta
    ( const DistMatrix<Int,U,GatheredDist<U>()>& perm,
      const DistMatrix<Int,U,GatheredDist<U>()>& invPerm );
};

// TODO: Generalize remaining routines to have prototypes of the form
//       (U,GatheredDist<U>())

// Apply column pivots
// ===================
template<typename T>
void ApplyColPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyColPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, 
  Int offset=0 );

template<typename T>
void ApplyInverseColPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyInverseColPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, 
  Int offset=0 );

// Apply row pivots
// ================
template<typename T>
void ApplyRowPivots( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyRowPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots,
  Int offset=0 );

template<typename T>
void ApplyInverseRowPivots
( Matrix<T>& A, const Matrix<Int>& pivots, Int offset=0 );
template<typename T,Dist U,Dist V,Dist UPerm>
void ApplyInverseRowPivots
( DistMatrix<T,U,V>& A, const DistMatrix<Int,UPerm,STAR>& pivots, 
  Int offset=0 );

// Apply symmetric pivots
// ======================
template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate=false, 
  Int offset=0 );
template<typename T,Dist UPerm>
void ApplySymmetricPivots
( UpperOrLower uplo, DistMatrix<T>& A,
  const DistMatrix<Int,UPerm,STAR>& pivots, bool conjugate=false, 
  Int offset=0 );

template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, Matrix<T>& A, 
  const Matrix<Int>& p, bool conjugate=false, 
  Int offset=0 );
template<typename T,Dist UPerm>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo, DistMatrix<T>& A,
  const DistMatrix<Int,UPerm,STAR>& pivots, bool conjugate=false, 
  Int offset=0 );

// Explicit permutation
// ====================
void ExplicitPermutation( const Matrix<Int>& perm, Matrix<Int>& P );
template<Dist UPerm,Dist U,Dist V>
void ExplicitPermutation
( const DistMatrix<Int,UPerm,STAR>& perm, DistMatrix<Int,U,V>& P );

// Invert permutation
// ==================
void InvertPermutation( const Matrix<Int>& perm, Matrix<Int>& invPerm );
template<Dist U>
void InvertPermutation
( const DistMatrix<Int,U,GatheredDist<U>()>& perm, 
        DistMatrix<Int,U,GatheredDist<U>()>& invPerm );

// Parity of a permutation
// =======================
bool PermutationParity( const Matrix<Int>& origPerm );
template<Dist UPerm>
bool PermutationParity( const DistMatrix<Int,UPerm,STAR>& origPerm );

// Permute columns
// ===============
template<typename T>
void PermuteCols( Matrix<T>& A, const Matrix<Int>& perm );
template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteCols
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm );

template<typename T>
void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& perm );
template<typename T,Dist U,Dist V,Dist UPerm>
void InversePermuteCols
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm );

template<typename T>
void PermuteCols
( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm );
template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteCols
( DistMatrix<T,U,V>& A,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm );

template<typename T,Dist U,Dist V>
void PermuteCols( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta );

// Permute rows
// ============
template<typename T>
void PermuteRows( Matrix<T>& A, const Matrix<Int>& perm );
template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteRows
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm );

template<typename T>
void InversePermuteRows( Matrix<T>& A, const Matrix<Int>& perm );
template<typename T,Dist U,Dist V,Dist UPerm>
void InversePermuteRows
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm );

template<typename T>
void PermuteRows
( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm );
template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteRows
( DistMatrix<T,U,V>& A,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm );

template<typename T,Dist U,Dist V>
void PermuteRows( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta );

// Parity of a sequence of pivots
// ==============================
bool PivotParity( const Matrix<Int>& p, Int pivotOffset=0 );
// TODO: Generalize implementation?
bool PivotParity( const DistMatrix<Int,VC,STAR>& p, Int pivotOffset=0 );

// Convert a pivot sequence to a partial permutation vector
// ========================================================
void PivotsToPartialPermutation
( const Matrix<Int>& pivots, Matrix<Int>& perm, Matrix<Int>& invPerm,
  Int offset=0 );
template<Dist U,Dist UPerm>
void PivotsToPartialPermutation
( const DistMatrix<Int,U,    STAR>& pivots,
        DistMatrix<Int,UPerm,STAR>& perm,
        DistMatrix<Int,UPerm,STAR>& invPerm, Int offset=0 );

// Convert a pivot sequence to a permutation vector
// ================================================
void PivotsToPermutation
( const Matrix<Int>& pivots, Matrix<Int>& perm, Int offset=0 );
template<Dist U,Dist UPerm>
void PivotsToPermutation
( const DistMatrix<Int,U,STAR>& pivots, DistMatrix<Int,UPerm,STAR>& perm,
  Int offset=0 );

void PivotsToInversePermutation
( const Matrix<Int>& pivots, Matrix<Int>& invPerm, Int offset=0 );
template<Dist U,Dist UPerm>
void PivotsToInversePermutation
( const DistMatrix<Int,U,STAR>& pivots, DistMatrix<Int,UPerm,STAR>& invPerm,
  Int offset=0 );

} // namespace El

#endif // ifndef EL_PERM_HPP
