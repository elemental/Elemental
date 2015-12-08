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

    parity_ = false;
    staleParity_ = false;

    swapSequence_ = true;

    // Only used if swapSequence_ = true
    // ---------------------------------
    nextSwapIndex_ = 0;
    implicitSwapOrigins_ = true;
    swapDests_.Empty();
    swapOrigins_.Empty();

    // Only used if swapSequence_ = false
    // ----------------------------------
    perm_.Empty();
    invPerm_.Empty();
    staleInverse_ = false;
}

void Permutation::ReserveSwaps( Int numSwaps )
{
    DEBUG_ONLY(CSE cse("Permutation::ReserveSwaps"))
    Empty();

    swapDests_.Resize( numSwaps, 1 );
    for( Int j=0; j<numSwaps; ++j )
        swapDests_.Set( j, 0, j );
}

void Permutation::AppendSwap( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("Permutation::AppendSwap"))
    if( !swapSequence_ )
        LogicError("Appending swaps to general permutations not yet enabled");
    if( !staleParity_ && origin != dest )
        parity_ = !parity_;

    if( !implicitSwapOrigins_ )
    {
        swapOrigins_.Set( nextSwapIndex_, 0, origin );
        swapDests_.Set( nextSwapIndex_, 0, dest );
        ++nextSwapIndex_;
        return;
    }

    if( origin == nextSwapIndex_ )
    {
        swapDests_.Set( nextSwapIndex_, 0, dest );
        ++nextSwapIndex_;
    }
    else
    {
        implicitSwapOrigins_ = false;
        const Int numSwaps = swapDests_.Height();

        swapOrigins_.Resize( numSwaps, 1 );
        for( Int j=0; j<numSwaps; ++j )
            swapOrigins_.Set( j, 0, j ); 

        swapOrigins_.Set( nextSwapIndex_, 0, origin );
        swapDests_.Set( nextSwapIndex_, 0, dest );
        ++nextSwapIndex_;
    }
}

void Permutation::AppendSwapSequence( const Permutation& perm, Int offset )
{
    DEBUG_ONLY(CSE cse("Permutation::AppendSwapSequence"))
    if( !swapSequence_ || !perm.swapSequence_ )
        LogicError("Appending swaps to general permutations not yet enabled");
    
    const Int numSwapAppends = perm.swapDests_.Height();
    if( perm.implicitSwapOrigins_ )
    {
        for( Int j=0; j<numSwapAppends; ++j )
            AppendSwap( nextSwapIndex_+j, perm.swapDests_.Get(j,0)-offset );
    }
    else
    {
        for( Int j=0; j<numSwapAppends; ++j )
            AppendSwap
            ( perm.swapOrigins_.Get(j,0)-offset,
              perm.swapDests_.Get(j,0)-offset );
    }
    nextSwapIndex_ += numSwapAppends;
}

void Permutation::SetImage( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("Permutation::SetImage"))
    if( swapSequence_ )
        LogicError("Cannot explicitly set the image of a swap sequence");
    perm_.Set( origin, 0, dest );
    invPerm_.Set( dest, 0, origin );
    staleParity_ = true;
}

void Permutation::MakeArbitrary( Int domainSize )
{
    DEBUG_ONLY(CSE cse("Permutation::MakeArbitrary"))
    if( swapSequence_ )
    {
        Empty();
        swapSequence_ = false;
    }
    perm_.Resize( domainSize, 1 );
    invPerm_.Resize( domainSize, 1 );
    for( Int j=0; j<domainSize; ++j )
        perm_.Set( j, 0, j );
    for( Int j=0; j<domainSize; ++j ) 
        invPerm_.Set( j, 0, j );
    parity_ = false;
    staleParity_ = false;
    staleInverse_ = false;
}

const Permutation& Permutation::operator=( const Permutation& p )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))

    swapSequence_ = p.swapSequence_;
    nextSwapIndex_ = p.nextSwapIndex_;
    implicitSwapOrigins_ = p.implicitSwapOrigins_;

    swapDests_ = p.swapDests_;
    swapOrigins_ = p.swapOrigins_;
    perm_ = p.perm_;
    invPerm_ = p.invPerm_;

    parity_ = p.parity_;
    staleParity_ = p.staleParity_;
    staleInverse_ = p.staleInverse_;

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

bool Permutation::IsSwapSequence() const
{ return swapSequence_; }

bool Permutation::IsImplicitSwapSequence() const
{ return implicitSwapOrigins_; }

template<typename T>
void Permutation::PermuteCols
( Matrix<T>& A, bool inverse, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteCols"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        const Int numSwaps = swapDests_.Height();
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }
            else
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }

        }
        else
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }
            else
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
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
        // TODO: Move El::PermuteCols into this class
        if( inverse )
            El::PermuteCols( A, invPerm_, perm_ );
        else
            El::PermuteCols( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::PermuteRows
( Matrix<T>& A, bool inverse, Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteRows"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        const Int numSwaps = swapDests_.Height();
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }
            else
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }

        }
        else
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }
            else
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
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
        if( inverse )
            El::PermuteRows( A, invPerm_, perm_ );
        else
            El::PermuteRows( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::PermuteSymmetrically
( UpperOrLower uplo,
  Matrix<T>& A,
  bool conjugate,
  bool inverse,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("Permutation::PermuteSymmetrically"))
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        const Int numSwaps = swapDests_.Height();
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }
            else
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }

        }
        else
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = j;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }
            else
            {
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = swapOrigins_.Get(j,0)-offset;
                    const Int dest = swapDests_.Get(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
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
        LogicError("General symmetric permutations are not yet supported");
    }
}

void Permutation::Explicit( Matrix<Int>& P ) const
{
    DEBUG_ONLY(CSE cse("Permutation::Explicit"))
    if( swapSequence_ )
    {
        Matrix<Int> perm;
        if( implicitSwapOrigins_ )
            PivotsToPermutation( swapDests_, perm );
        else
            LogicError("Unsupported explicit permutation option");
        ExplicitPermutation( perm, P );
    }
    else
    {
        ExplicitPermutation( perm_, P );
    }
}

#define PROTO(T) \
  template void Permutation::PermuteCols \
  ( Matrix<T>& A, \
    bool inverse, \
    Int offset ) const; \
  template void Permutation::PermuteRows \
  ( Matrix<T>& A, \
    bool inverse, \
    Int offset ) const; \
  template void Permutation::PermuteSymmetrically \
  ( UpperOrLower uplo, \
    Matrix<T>& A, \
    bool conjugate, \
    bool inverse, \
    Int offset ) const;

#include "El/macros/Instantiate.h"

} // namespace El
