/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PERM_PERMUTATION_HPP
#define EL_PERM_PERMUTATION_HPP

namespace El {

class Permutation
{
public:
    Permutation();

    void Empty();
    void MakeIdentity( Int size );

    void ReserveSwaps( Int maxSwaps );
    void MakeArbitrary() const;

    const Permutation& operator=( const Permutation& P );

    void Swap( Int origin, Int dest );
    void SwapSequence( const Permutation& P, Int offset=0 );
    void SwapSequence
    ( const Matrix<Int>& swapOrigins,
      const Matrix<Int>& swapDests, Int offset=0 );

    void ImplicitSwapSequence( const Matrix<Int>& swapDests, Int offset=0 );

    // Explicit image queries or modifications force the permutation to switch
    // from an explicit sequence of swaps to storing an explicit permutation
    Int Image( Int origin ) const;
    Int Preimage( Int dest ) const;
    void SetImage( Int origin, Int dest );

    // The following return the same result but follow the usual convention
    Int Height() const;
    Int Width() const;

    bool Parity() const;
    bool IsSwapSequence() const;
    bool IsImplicitSwapSequence() const;

    // NOTE: This is only valid if IsImplicitSwapSequence() is true, otherwise
    //       it is implicitly [0,...,numSwaps-1]
    const Matrix<Int> SwapOrigins() const;
    // NOTE: This is only valid if IsSwapSequence() is true
    const Matrix<Int> SwapDestinations() const;

    template<typename T>
    void PermuteCols
    ( Matrix<T>& A,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteCols
    ( Matrix<T>& A,
      Int offset=0 ) const;

    template<typename T>
    void PermuteRows
    ( Matrix<T>& A,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteRows
    ( Matrix<T>& A,
      Int offset=0 ) const;

    template<typename T>
    void PermuteSymmetrically
    ( UpperOrLower uplo,
      Matrix<T>& A,
      bool conjugate=false,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteSymmetrically
    ( UpperOrLower uplo,
      Matrix<T>& A,
      bool conjugate=false,
      Int offset=0 ) const;

    // Form the permutation vector p so that P A = A(p,:)
    void ExplicitVector( Matrix<Int>& p ) const;

    // Form the permutation matrix P so that P A = A(p,:)
    void ExplicitMatrix( Matrix<Int>& P ) const;

private:

    Int size_=0;

    mutable bool parity_=false;
    mutable bool staleParity_=false;

    mutable bool swapSequence_=true;

    // Only used if swapSequence_=true
    // -------------------------------
    // NOTE: As swaps are added, if the origin sequence is of the usual form
    //       of 0, 1, 2, ..., then an explicit swap origin vector is not
    //       maintained. However, if an unexpected origin is ever encountered,
    //       then an explicit list is then maintained.
    mutable Int numSwaps_=0;
    mutable bool implicitSwapOrigins_=true;
    mutable Matrix<Int> swapDests_, swapOrigins_;

    // Only used if swapSequence_=false
    // --------------------------------
    mutable Matrix<Int> perm_;
    mutable Matrix<Int> invPerm_;
    mutable bool staleInverse_=true;

    friend class DistPermutation;
};

} // namespace El

#endif // ifndef EL_PERM_PERMUTATION_HPP
