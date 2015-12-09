/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

Permutation::Permutation() { }

void Permutation::Empty()
{
    DEBUG_ONLY(CSE cse("Permutation::Empty"))
    size_ = 0;

    parity_ = false;
    staleParity_ = false;

    swapSequence_ = true;

    // Only used if swapSequence_ = true
    // ---------------------------------
    numSwaps_ = 0;
    implicitSwapOrigins_ = true;
    swapDests_.Empty();
    swapOrigins_.Empty();

    // Only used if swapSequence_ = false
    // ----------------------------------
    perm_.Empty();
    invPerm_.Empty();
    staleInverse_ = false;
}

void Permutation::MakeIdentity( Int size )
{
    DEBUG_ONLY(CSE cse("Permutation::MakeIdentity"))

    if( !swapSequence_ )
    {
        Empty();
        size_ = size;
        return;
    }

    // Carefully ensure that a previous swap reservation is not blocked
    // if the permutation is still in swap mode

    size_ = size;

    parity_ = false;
    staleParity_ = false;

    numSwaps_ = 0;
    implicitSwapOrigins_ = true;
}

void Permutation::ReserveSwaps( Int maxSwaps )
{
    DEBUG_ONLY(CSE cse("Permutation::ReserveSwaps"))

    // Arbitrary permutations can trivially support arbitrarily many swaps
    if( !swapSequence_ )
        return;

    if( maxSwaps > swapDests_.Height() )
    {
        MakeIdentity( size_ );
        swapDests_.Resize( maxSwaps, 1 );
        return;
    }
}

void Permutation::RowSwap( Int origin, Int dest )
{
    DEBUG_ONLY(
      CSE cse("Permutation::RowSwap");
      if( origin < 0 || origin >= size_ || dest < 0 || dest >= size_ )
          LogicError
          ("Attempted swap (",origin,",",dest,") for perm. of size ",size_);
    )
    if( origin != dest )
        parity_ = !parity_;
    if( swapSequence_ && numSwaps_ == swapDests_.Height() )
        MakeArbitrary();
    if( !swapSequence_ )
    {
        El::RowSwap( perm_, origin, dest );
        staleInverse_ = true;
        return;
    }

    if( !implicitSwapOrigins_ )
    {
        swapOrigins_.Set( numSwaps_, 0, origin );
        swapDests_.Set( numSwaps_, 0, dest );
        ++numSwaps_;
        return;
    }

    if( origin == numSwaps_ )
    {
        swapDests_.Set( numSwaps_, 0, dest );
        ++numSwaps_;
    }
    else
    {
        implicitSwapOrigins_ = false;
        const Int numSwaps = swapDests_.Height();

        swapOrigins_.Resize( numSwaps, 1 );
        for( Int j=0; j<numSwaps; ++j )
            swapOrigins_.Set( j, 0, j ); 

        swapOrigins_.Set( numSwaps_, 0, origin );
        swapDests_.Set( numSwaps_, 0, dest );
        ++numSwaps_;
    }
}

void Permutation::RowSwapSequence( const Permutation& P, Int offset )
{
    DEBUG_ONLY(CSE cse("Permutation::RowSwapSequence"))
    if( P.swapSequence_ )
    {
        const Int numSwapAppends = P.numSwaps_;
        if( P.implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwapAppends; ++j )
                RowSwap( j+offset, P.swapDests_.Get(j,0)+offset );
        }
        else
        {
            for( Int j=0; j<numSwapAppends; ++j )
                RowSwap
                ( P.swapOrigins_.Get(j,0)+offset,
                  P.swapDests_.Get(j,0)+offset );
        }
    }
    else
    {
        MakeArbitrary();
        P.PermuteRows( perm_, offset );
        invPerm_.Empty();

        staleParity_ = true;
        staleInverse_ = true;
    }
}

void Permutation::SetImage( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("Permutation::SetImage"))
    MakeArbitrary();
    perm_.Set( origin, 0, dest );
    staleInverse_ = true;
    staleParity_ = true;
}

void Permutation::MakeArbitrary()
{
    DEBUG_ONLY(CSE cse("Permutation::MakeArbitrary"))
    if( !swapSequence_ )
        return;

    // Compose the swaps into an explicit permutation vector
    // -----------------------------------------------------
    perm_.Resize( size_, 1 );
    for( Int j=0; j<size_; ++j )
        perm_.Set( j, 0, j );
    PermuteRows( perm_ );

    invPerm_.Empty();
    staleInverse_ = true;

    // Clear the swap information
    // --------------------------
    swapSequence_ = false;
    numSwaps_ = 0;
    implicitSwapOrigins_ = true;
    swapDests_.Empty();
    swapOrigins_.Empty();
}

const Permutation& Permutation::operator=( const Permutation& P )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))

    size_ = P.size_;

    swapSequence_ = P.swapSequence_;
    numSwaps_ = P.numSwaps_;
    implicitSwapOrigins_ = P.implicitSwapOrigins_;

    swapDests_ = P.swapDests_;
    swapOrigins_ = P.swapOrigins_;
    perm_ = P.perm_;
    invPerm_ = P.invPerm_;

    parity_ = P.parity_;
    staleParity_ = P.staleParity_;
    staleInverse_ = P.staleInverse_;

    return *this;
}

bool Permutation::Parity() const
{
    DEBUG_ONLY(CSE cse("Permutation::Parity"))
    if( staleParity_ )
    {
        LogicError("Computation of parity not yet enabled");
        staleParity_ = false;
        return parity_; 
    }
    return parity_;
}

Int Permutation::Height() const
{ return size_; }

Int Permutation::Width() const
{ return size_; }

bool Permutation::IsSwapSequence() const
{ return swapSequence_; }

bool Permutation::IsImplicitSwapSequence() const
{ return implicitSwapOrigins_; }

template<typename T>
void Permutation::PermuteCols( Matrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteCols"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::InversePermuteCols( Matrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::InversePermuteCols"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, invPerm_, perm_ );
    }
}

template<typename T>
void Permutation::PermuteRows( Matrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteRows"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
    }
    else
    {
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteRows into this class
        El::PermuteRows( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::InversePermuteRows( Matrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::InversePermuteRows"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteRows into this class
        El::PermuteRows( A, invPerm_, perm_ );
    }
}

template<typename T>
void Permutation::PermuteSymmetrically
( UpperOrLower uplo,
  Matrix<T>& A,
  bool conjugate,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteSymmetrically"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        LogicError("General symmetric permutations are not yet supported");
    }
}

template<typename T>
void Permutation::InversePermuteSymmetrically
( UpperOrLower uplo,
  Matrix<T>& A,
  bool conjugate,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::InversePermuteSymmetrically"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_.Get(j,0)+offset;
                const Int dest = swapDests_.Get(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        LogicError("General symmetric permutations are not yet supported");
    }
}

void Permutation::Explicit( Matrix<Int>& P ) const
{
    DEBUG_ONLY(CSE cse("Permutation::Explicit"))
    // TODO: Faster algorithm
    Identity( P, size_, size_ );
    PermuteRows( P );
}

#define PROTO(T) \
  template void Permutation::PermuteCols \
  ( Matrix<T>& A, \
    Int offset ) const; \
  template void Permutation::InversePermuteCols \
  ( Matrix<T>& A, \
    Int offset ) const; \
  template void Permutation::PermuteRows \
  ( Matrix<T>& A, \
    Int offset ) const; \
  template void Permutation::InversePermuteRows \
  ( Matrix<T>& A, \
    Int offset ) const; \
  template void Permutation::PermuteSymmetrically \
  ( UpperOrLower uplo, \
    Matrix<T>& A, \
    bool conjugate, \
    Int offset ) const; \
  template void Permutation::InversePermuteSymmetrically \
  ( UpperOrLower uplo, \
    Matrix<T>& A, \
    bool conjugate, \
    Int offset ) const;

#include "El/macros/Instantiate.h"

} // namespace El
