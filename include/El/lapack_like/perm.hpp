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

// TODO: Maintain a hash map from row/column alignments to instances of this
//       class within DistPermutation for both row and column permutations
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
    ( const DistMatrix<Int,STAR,STAR>& p,
      const DistMatrix<Int,STAR,STAR>& pInv,
            Int permAlign,
            mpi::Comm permComm );

    void Update
    ( const DistMatrix<Int,STAR,STAR>& p,
      const DistMatrix<Int,STAR,STAR>& pInv,
            Int permAlign,
            mpi::Comm permComm );
};

class Permutation
{
public:
    Permutation();

    void Resize( Int domainSize );
    void AppendSwap( Int origin, Int dest );
    void AppendSwapSequence( const Permutation& perm, Int offset=0 );

    bool Parity() const;
    bool IsSwapSequence() const;
    bool IsImplicitSwapSequence() const;

    void Explicit( Matrix<Int>& P ) const;

    template<typename T>
    void ApplyToCols
    ( Matrix<T>& A,
      bool inverse=false,
      Int offset=0 );

    template<typename T>
    void ApplyToRows
    ( Matrix<T>& A,
      bool inverse=false,
      Int offset=0 );

    template<typename T>
    void ApplySymmetrically
    ( UpperOrLower uplo,
      Matrix<T>& A,
      bool conjugate=false,
      bool inverse=false,
      Int offset=0 );

private:
    Int domainSize_=0;
    mutable bool parity_=false;
    mutable bool staleParity_=false;

    bool swapSequence_=true;

    // Only used if swapSequence=true
    // ------------------------------
    // NOTE: As swaps are added, if the origin sequence is of the usual form
    //       of 0, 1, 2, ..., then an explicit swap origin vector is not
    //       maintained. However, if an unexpected origin is ever encountered,
    //       then an explicit list is then maintained.
    Int swapOffset_=0;
    Int lastSwapIndex_=-1;
    bool implicitSwapOrigins_=true;
    Matrix<Int> swapDests_, swapOrigins_;

    // Only used if swapSequence=false
    // -------------------------------
            Matrix<Int> perm_;
    mutable Matrix<Int> invPerm_;
    mutable bool staleInverse_=false;
};

class DistPermutation
{
public:
    DistPermutation( const Grid& g );

    void Resize( Int domainSize );
    void AppendSwap( Int origin, Int dest );
    void AppendSwapSequence( const DistPermutation& perm, Int offset=0 );

    bool Parity() const;
    bool IsSwapSequence() const;
    bool IsImplicitSwapSequence() const;

    void Explicit( AbstractDistMatrix<Int>& P ) const;

    template<typename T>
    void ApplyToCols
    ( AbstractDistMatrix<T>& A,
      bool inverse=false,
      Int offset=0 );

    template<typename T>
    void ApplyToRows
    ( AbstractDistMatrix<T>& A,
      bool inverse=false,
      Int offset=0 );

    template<typename T>
    void ApplySymmetrically
    ( UpperOrLower uplo,
      AbstractDistMatrix<T>& A,
      bool conjugate=false,
      bool inverse=false,
      Int offset=0 );

    //void FormMetadata() const;

private:
    Int domainSize_=0;
    mutable bool parity_=false;
    mutable bool staleParity_=false;

    //bool staleMetadata_=false;
    //mutable PermutationMeta meta_;

    bool swapSequence_=true;

    // Only used if swapSequence=true
    // ------------------------------
    // NOTE: As swaps are added, if the origin sequence is of the usual form
    //       of 0, 1, 2, ..., then an explicit swap origin vector is not
    //       maintained. However, if an unexpected origin is ever encountered,
    //       then an explicit list is then maintained.
    Int swapOffset_=0;
    Int lastSwapIndex_=-1;
    bool implicitSwapOrigins_=true;
    DistMatrix<Int,VC,STAR> swapDests_, swapOrigins_;

    // Only used if swapSequence=false
    // -------------------------------
            DistMatrix<Int,VC,STAR> perm_;
    mutable DistMatrix<Int,VC,STAR> invPerm_;
    mutable bool staleInverse_=false;
};

// Permute columns
// ===============
// TODO: Deprecate this
template<typename T>
void PermuteCols
(       AbstractDistMatrix<T>& A,
  const PermutationMeta& meta,
  bool inverse=false ); 

// Permute rows
// ============
// TODO: Deprecate this
template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const PermutationMeta& meta,
  bool inverse=false ); 

// Symmetric permutation
// =====================
template<typename T>
void SymmetricPermutation
( UpperOrLower uplo,
        AbstractDistMatrix<T>& A,
  const DistPermutation& perm,
  bool inverse=false,
  bool conjugate=false );

// Deprecation line
// ################

// Apply column pivots
// ===================
// NOTE: These routine are now deprecated

template<typename T>
void ApplyColPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset=0 );
template<typename T>
void ApplyColPivots
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots, 
        Int offset=0 );

template<typename T>
void ApplyInverseColPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset=0 );
template<typename T>
void ApplyInverseColPivots
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots, 
        Int offset=0 );

// Apply row pivots
// ================
// NOTE: These routine are now deprecated

template<typename T>
void ApplyRowPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset=0 );
template<typename T>
void ApplyRowPivots
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots,
        Int offset=0 );

template<typename T>
void ApplyInverseRowPivots
(       Matrix<T>& A,
  const Matrix<Int>& pivots,
        Int offset=0 );
template<typename T>
void ApplyInverseRowPivots
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots, 
        Int offset=0 );

// Apply symmetric pivots
// ======================
// NOTE: These routine are now deprecated
template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo,
        Matrix<T>& A, 
  const Matrix<Int>& p,
  bool conjugate=false, 
  Int offset=0 );
template<typename T>
void ApplySymmetricPivots
( UpperOrLower uplo,
        AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots,
  bool conjugate=false, 
  Int offset=0 );

template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo,
        Matrix<T>& A, 
  const Matrix<Int>& p,
  bool conjugate=false, 
  Int offset=0 );
template<typename T>
void ApplyInverseSymmetricPivots
( UpperOrLower uplo,
        AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& pivots,
  bool conjugate=false, 
  Int offset=0 );

// Explicit permutation
// ====================
// NOTE: These routine are now deprecated
void ExplicitPermutation
( const Matrix<Int>& p,
        Matrix<Int>& P );
void ExplicitPermutation
( const AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& P );

// Invert permutation
// ==================
// NOTE: These routine are now deprecated
void InvertPermutation
( const Matrix<Int>& p,
        Matrix<Int>& pInv );
void InvertPermutation
( const AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv );

// Parity of a permutation
// =======================
// NOTE: These routine are now deprecated
bool PermutationParity( const Matrix<Int>& p );
bool PermutationParity( const AbstractDistMatrix<Int>& p );

// Permute columns
// ===============
// NOTE: These routine are now deprecated
template<typename T>
void PermuteCols
(       Matrix<T>& A,
  const Matrix<Int>& p );
template<typename T>
void PermuteCols
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p );

template<typename T>
void InversePermuteCols
(       Matrix<T>& A,
  const Matrix<Int>& p );
template<typename T>
void InversePermuteCols
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p );

template<typename T>
void PermuteCols
(       Matrix<T>& A,
  const Matrix<Int>& p,
  const Matrix<Int>& pInv );
template<typename T>
void PermuteCols
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p,
  const AbstractDistMatrix<Int>& pInv );

// Permute rows
// ============
// NOTE: These routine are now deprecated
template<typename T>
void PermuteRows
(       Matrix<T>& A,
  const Matrix<Int>& p );
template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p );

template<typename T>
void InversePermuteRows
(       Matrix<T>& A, 
  const Matrix<Int>& p );
template<typename T>
void InversePermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p );

template<typename T>
void PermuteRows
(       Matrix<T>& A,
  const Matrix<Int>& p,
  const Matrix<Int>& pInv );
template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& p,
  const AbstractDistMatrix<Int>& pInv );

// Parity of a sequence of pivots
// ==============================
// NOTE: These routine are now deprecated
bool PivotParity( const Matrix<Int>& p, Int pivotOffset=0 );
bool PivotParity( const AbstractDistMatrix<Int>& p, Int pivotOffset=0 );

// Convert a pivot sequence to a partial permutation vector
// ========================================================
// NOTE: These routine are now deprecated
void PivotsToPartialPermutation
( const Matrix<Int>& pivots,
        Matrix<Int>& p,
        Matrix<Int>& pInv,
        Int offset=0 );
void PivotsToPartialPermutation
( const AbstractDistMatrix<Int>& pivots,
        AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv,
        Int offset=0 );
void PivotsToPartialPermutation
( const DistMatrix<Int,STAR,STAR>& pivots,
        AbstractDistMatrix<Int>& p,
        AbstractDistMatrix<Int>& pInv,
        Int offset=0 );

// Convert a pivot sequence to a permutation vector
// ================================================
// NOTE: These routine are now deprecated
void PivotsToPermutation
( const Matrix<Int>& pivots,
        Matrix<Int>& p,
        Int offset=0 );
void PivotsToPermutation
( const AbstractDistMatrix<Int>& pivots,
        AbstractDistMatrix<Int>& p,
        Int offset=0 );

void PivotsToInversePermutation
( const Matrix<Int>& pivots,
        Matrix<Int>& pInv,
        Int offset=0 );
void PivotsToInversePermutation
( const AbstractDistMatrix<Int>& pivots,
        AbstractDistMatrix<Int>& pInv,
        Int offset=0 );

} // namespace El

#endif // ifndef EL_PERM_HPP
