/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PERM_PERMUTATION_HPP
#define EL_PERM_PERMUTATION_HPP

namespace El {

class Permutation
{
public:
    Permutation();

    void Empty();
    void ReserveSwaps( Int numSwaps );
    void MakeArbitrary( Int domainSize );

    const Permutation& operator=( const Permutation& p );

    void AppendSwap( Int origin, Int dest );
    void AppendSwapSequence( const Permutation& perm, Int offset=0 );

    void SetImage( Int origin, Int dest );

    bool Parity() const;
    bool IsSwapSequence() const;
    bool IsImplicitSwapSequence() const;

    template<typename T>
    void PermuteCols
    ( Matrix<T>& A,
      bool inverse=false,
      Int offset=0 ) const;

    template<typename T>
    void PermuteRows
    ( Matrix<T>& A,
      bool inverse=false,
      Int offset=0 ) const;

    template<typename T>
    void PermuteSymmetrically
    ( UpperOrLower uplo,
      Matrix<T>& A,
      bool conjugate=false,
      bool inverse=false,
      Int offset=0 ) const;

    void Explicit( Matrix<Int>& P ) const;

private:

    mutable bool parity_=false;
    mutable bool staleParity_=false;

    bool swapSequence_=true;

    // Only used if swapSequence_=true
    // -------------------------------
    // NOTE: As swaps are added, if the origin sequence is of the usual form
    //       of 0, 1, 2, ..., then an explicit swap origin vector is not
    //       maintained. However, if an unexpected origin is ever encountered,
    //       then an explicit list is then maintained.
    Int nextSwapIndex_=0;
    bool implicitSwapOrigins_=true;
    Matrix<Int> swapDests_, swapOrigins_;

    // Only used if swapSequence_=false
    // --------------------------------
            Matrix<Int> perm_;
    mutable Matrix<Int> invPerm_;
    mutable bool staleInverse_=true;

    friend class DistPermutation;
};

} // namespace El

#endif // ifndef EL_PERM_PERMUTATION_HPP
