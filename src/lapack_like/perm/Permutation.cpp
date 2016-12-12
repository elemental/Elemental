/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace {

template<typename T>
void PermuteCols
(       Matrix<T>& A,
  const Matrix<Int>& perm,
  const Matrix<Int>& invPerm )
{
    EL_DEBUG_CSE
    const Int b = perm.Height();
    EL_DEBUG_ONLY(
      if( A.Width() < b || b != invPerm.Height() )
          LogicError
          ("perm and invPerm must be vectors of equal length that are not "
           "wider than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // Make a copy of the first b columns
    auto AColPanView = A( IR(0,m), IR(0,b) );
    auto AColPanCopy( AColPanView );

    // Make a copy of the preimage columns
    Matrix<T> APreimageCopy( m, b );

    for( Int j=0; j<b; ++j )
    {
        const Int jPre = perm(j);
        if( jPre >= b )
            MemCopy( &APreimageCopy(0,j), &A(0,jPre), m );
    }

    // Apply the permutations
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = perm(j);
        const Int jPost = invPerm(j);
        // Move row[j] into row[jPost]
        MemCopy( &A(0,jPost), &AColPanCopy(0,j), m );
        // Move row[jPre] into row[j]
        if( jPre >= b )
            MemCopy( &A(0,j), &APreimageCopy(0,j), m );
    }
}

template<typename T>
void PermuteRows
(       Matrix<T>& A,
  const Matrix<Int>& perm,
  const Matrix<Int>& invPerm )
{
    EL_DEBUG_CSE
    const Int b = perm.Height();
    EL_DEBUG_ONLY(
      if( A.Height() < b || b != invPerm.Height() )
          LogicError
          ("perm and invPerm must be vectors of equal length that are not "
           "taller than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // Make a copy of the first b rows
    auto ARowPanView = A( IR(0,b), IR(0,n) );
    Matrix<T> ARowPanTrans;
    Transpose( ARowPanView, ARowPanTrans );

    // Make a copy of the preimage rows
    Matrix<T> APreimageTrans( n, b );

    for( Int i=0; i<b; ++i )
    {
        const Int iPre = perm(i);
        if( iPre >= b )
            for( Int j=0; j<n; ++j )
                APreimageTrans(j,i) = A(iPre,j);
    }

    // Apply the permutations
    for( Int i=0; i<b; ++i )
    {
        const Int iPre = perm(i);
        const Int iPost = invPerm(i);
        // Move row[i] into row[image[i]]
        for( Int j=0; j<n; ++j )
            A(iPost,j) = ARowPanTrans(j,i);
        if( iPre >= b )
        {
            // Move row[preimage[i]] into row[i]
            for( Int j=0; j<n; ++j )
                A(i,j) = APreimageTrans(j,i);
        }
    }
}

void InvertPermutation
( const Matrix<Int>& p, Matrix<Int>& pInv )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( p.Width() != 1 )
          LogicError("p must be a column vector");
    )
    const Int n = p.Height();
    pInv.Resize( n, 1 );
    if( n == 0 )
        return;

    EL_DEBUG_ONLY(
      // This is obviously necessary but not sufficient for 'p' to contain
      // a reordering of (0,1,...,n-1).
      const Int range = MaxNorm( p ) + 1;
      if( range != n )
          LogicError("Invalid permutation range");
    )

    for( Int i=0; i<n; ++i )
        pInv(p(i)) = i;
}

} // anonymous namespace

Permutation::Permutation() { }

void Permutation::Empty()
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE

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
    EL_DEBUG_CSE

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

void Permutation::Swap( Int origin, Int dest )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
        swapOrigins_(numSwaps_) = origin;
        swapDests_(numSwaps_) = dest;
        ++numSwaps_;
        return;
    }

    if( origin == numSwaps_ )
    {
        swapDests_(numSwaps_) = dest;
        ++numSwaps_;
    }
    else
    {
        implicitSwapOrigins_ = false;
        const Int maxSwaps = swapDests_.Height();

        swapOrigins_.Resize( maxSwaps, 1 );
        for( Int j=0; j<maxSwaps; ++j )
            swapOrigins_(j) = j;

        swapOrigins_(numSwaps_) = origin;
        swapDests_(numSwaps_) = dest;
        ++numSwaps_;
    }
}

void Permutation::SwapSequence( const Permutation& P, Int offset )
{
    EL_DEBUG_CSE
    if( P.swapSequence_ )
    {
        const Int numSwapAppends = P.numSwaps_;
        if( P.implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwapAppends; ++j )
                Swap( j+offset, P.swapDests_(j)+offset );
        }
        else
        {
            for( Int j=0; j<numSwapAppends; ++j )
                Swap( P.swapOrigins_(j)+offset, P.swapDests_(j)+offset );
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

void Permutation::SwapSequence
( const Matrix<Int>& swapOrigins,
  const Matrix<Int>& swapDests,
  Int offset )
{
    EL_DEBUG_CSE
    // TODO: Assert swapOrigins and swapDests are column vectors of same size
    const Int numSwaps = swapDests.Height();
    for( Int k=0; k<numSwaps; ++k )
        Swap( swapOrigins(k)+offset, swapDests(k)+offset );
}

void Permutation::ImplicitSwapSequence
( const Matrix<Int>& swapDests,
  Int offset )
{
    EL_DEBUG_CSE
    const Int numPrevSwaps = numSwaps_;

    // TODO: Assert swapOrigins and swapDests are column vectors of same size
    const Int numSwaps = swapDests.Height();
    for( Int k=0; k<numSwaps; ++k )
        Swap( numPrevSwaps+k, swapDests(k)+offset );
}

Int Permutation::Image( Int origin ) const
{
    EL_DEBUG_CSE
    MakeArbitrary();
    if( staleInverse_ )
    {
        El::InvertPermutation( perm_, invPerm_ );
        staleInverse_ = false;
    }
    return invPerm_(origin);
}

Int Permutation::Preimage( Int dest ) const
{
    EL_DEBUG_CSE
    MakeArbitrary();
    return perm_(dest);
}

void Permutation::SetImage( Int origin, Int dest )
{
    EL_DEBUG_CSE
    MakeArbitrary();
    perm_(dest) = origin;
    invPerm_(origin) = dest;
    staleParity_ = true;
}

void Permutation::MakeArbitrary() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ )
        return;

    // Compose the swaps into an explicit permutation vector
    // -----------------------------------------------------
    perm_.Resize( size_, 1 );
    for( Int j=0; j<size_; ++j )
        perm_(j) = j;
    PermuteRows( perm_ );

    invPerm_.Resize( size_, 1 );
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
    EL_DEBUG_CSE

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
    EL_DEBUG_CSE
    if( staleParity_ )
    {
        if( swapSequence_ )
            LogicError("Unexpected stale parity for a swap sequence");

        // Walking through the process of LU factorization with partial 
        // pivoting for a permutation matrix, which never requires a
        // Schur-complement update, yields an algorithm for expressing the 
        // inverse of a permutation in terms of a sequence of transpositions in         // linear time. Note that performing the swaps requires access to the 
        // inverse permutation, which can be formed in linear time.
        if( staleInverse_ )
        {
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        auto permCopy( perm_ );
        auto invPermCopy( invPerm_ );

        parity_ = false;
        for( Int k=0; k<size_; ++k )
        {
            const Int permVal = permCopy(k);
            if( permVal != k )
            {
                parity_ = !parity_;
                const Int invPermVal = invPermCopy(k);
                // We only need to perform half of the swaps
                //      perm[k] <-> perm[invPerm[k]]
                //   invPerm[k] <-> invPerm[perk[k]] 
                // since we will not need to access perm[k] and invPerm[k] again
                permCopy(invPermVal) = permVal;
                invPermCopy(permVal) = invPermVal;
            }
        }

        staleParity_ = false;
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

const Matrix<Int> Permutation::SwapOrigins() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ || !implicitSwapOrigins_ )
        LogicError("Swap origins are not explicitly stored");
    return swapOrigins_(IR(0,numSwaps_),ALL);
}

const Matrix<Int> Permutation::SwapDestinations() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ )
        LogicError("Swap destinations are not explicitly stored");
    return swapDests_(IR(0,numSwaps_),ALL);
}

template<typename T>
void Permutation::PermuteCols( Matrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
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
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::InversePermuteCols( Matrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
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
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, invPerm_, perm_ );
    }
}

template<typename T>
void Permutation::PermuteRows( Matrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
    }
    else
    {
        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        // TODO: Move El::PermuteRows into this class
        El::PermuteRows( A, perm_, invPerm_ );
    }
}

template<typename T>
void Permutation::InversePermuteRows( Matrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
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
            El::InvertPermutation( perm_, invPerm_ );
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
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
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
            El::InvertPermutation( perm_, invPerm_ );
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
    EL_DEBUG_CSE
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
                const Int dest = swapDests_(j)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = swapOrigins_(j)+offset;
                const Int dest = swapDests_(j)+offset;
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
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        LogicError("General symmetric permutations are not yet supported");
    }
}

void Permutation::ExplicitVector( Matrix<Int>& p ) const
{
    EL_DEBUG_CSE
    if( swapSequence_ )
    {   
        p.Resize( size_, 1 );
        for( Int i=0; i<size_; ++i )
            p(i) = i;
        PermuteRows( p );
    }
    else
    {
        p = perm_;
    }
}

void Permutation::ExplicitMatrix( Matrix<Int>& P ) const
{
    EL_DEBUG_CSE
    Matrix<Int> p;
    ExplicitVector( p );

    Zeros( P, size_, size_ );
    for( Int i=0; i<size_; ++i )
        P(i,p(i)) = 1;
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

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
