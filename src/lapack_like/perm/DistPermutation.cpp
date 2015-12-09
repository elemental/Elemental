/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

DistPermutation::DistPermutation( const Grid& g )
: grid_(g)
{ 
    DEBUG_ONLY(CSE cse("DistPermutation::DistPermutation"))
    swapDests_.SetGrid( g );
    swapOrigins_.SetGrid( g );
    perm_.SetGrid( g );
    invPerm_.SetGrid( g );
}

void DistPermutation::SetGrid( const Grid& g )
{
    DEBUG_ONLY(CSE cse("DistPermutation::SetGrid"))
    Empty();
    swapDests_.SetGrid( g );
    swapOrigins_.SetGrid( g );
    perm_.SetGrid( g );
    invPerm_.SetGrid( g );
}

void DistPermutation::Empty()
{
    DEBUG_ONLY(CSE cse("DistPermutation::Empty"))
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

    rowMeta_.clear();
    colMeta_.clear();
    staleMeta_ = false;
}

void DistPermutation::MakeIdentity( Int size )
{
    DEBUG_ONLY(CSE cse("DistPermutation::MakeIdentity"))

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

void DistPermutation::ReserveSwaps( Int maxSwaps )
{
    DEBUG_ONLY(CSE cse("DistPermutation::ReserveSwaps"))

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

void DistPermutation::RowSwap( Int origin, Int dest )
{
    DEBUG_ONLY(
      CSE cse("DistPermutation::RowSwap");
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
        staleMeta_ = true;
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
        const Int maxSwaps = swapDests_.Height();

        swapOrigins_.SetGrid( grid_ );
        swapOrigins_.Resize( maxSwaps, 1 );
        
        auto swapOriginsActive = swapOrigins_(IR(0,numSwaps_),ALL);
        const Int localHeight = swapOriginsActive.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            swapOriginsActive.SetLocal
            ( iLoc, 0, swapOriginsActive.GlobalRow(iLoc) );

        swapOrigins_.Set( numSwaps_, 0, origin );
        swapDests_.Set( numSwaps_, 0, dest );
        ++numSwaps_;
    }
}

void DistPermutation::RowSwapSequence
( const DistPermutation& P, Int offset )
{
    DEBUG_ONLY(CSE cse("DistPermutation::RowSwapSequence"))
    if( P.swapSequence_ )
    {
        const Int numSwapAppends = P.numSwaps_;
        auto activeInd = IR(0,numSwapAppends);

        DistMatrix<Int,STAR,STAR> swapDests_STAR_STAR =
          P.swapDests_(activeInd,ALL);
        if( P.implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwapAppends; ++j )
                RowSwap
                ( j+offset, swapDests_STAR_STAR.GetLocal(j,0)+offset );
        }
        else
        {
            DistMatrix<Int,STAR,STAR> swapOrigins_STAR_STAR =
              P.swapOrigins_(activeInd,ALL);
            for( Int j=0; j<numSwapAppends; ++j )
                RowSwap
                ( swapOrigins_STAR_STAR.GetLocal(j,0)+offset,
                  swapDests_STAR_STAR.GetLocal(j,0)+offset );
        }
    }
    else
    {
        MakeArbitrary();
        P.PermuteRows( perm_, offset );
        invPerm_.Empty();

        staleParity_ = true;
        staleInverse_ = true;
        staleMeta_ = true; 
    }
}

void DistPermutation::SetImage( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("DistPermutation::SetImage"))
    MakeArbitrary();
    perm_.Set( origin, 0, dest );
    staleParity_ = true;
    staleInverse_ = true;
    staleMeta_ = true;
}

void DistPermutation::MakeArbitrary()
{
    DEBUG_ONLY(CSE cse("DistPermutation::MakeArbitrary"))
    if( !swapSequence_ )
        return;

    // Compose the swaps into an explicit permutation vector
    // -----------------------------------------------------
    perm_.Resize( size_, 1 );
    const Int localHeight = perm_.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        perm_.SetLocal( iLoc, 0, perm_.GlobalRow(iLoc) );
    PermuteRows( perm_ );

    invPerm_.Empty();
    staleInverse_ = true;
    staleMeta_ = true;

    // Clear the swap information
    // --------------------------
    swapSequence_ = false;
    numSwaps_ = 0;
    implicitSwapOrigins_ = true;
    swapDests_.Empty();
    swapOrigins_.Empty();
}

const DistPermutation& DistPermutation::operator=( const Permutation& P )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))
    if( grid_.Size() != 1 )
        LogicError("Invalid grid size for sequential copy");

    size_ = P.size_;

    swapSequence_ = P.swapSequence_;
    numSwaps_ = P.numSwaps_;
    implicitSwapOrigins_ = P.implicitSwapOrigins_;

    swapDests_.Resize( P.swapDests_.Height(), 1 );
    swapOrigins_.Resize( P.swapOrigins_.Height(), 1 );
    perm_.Resize( P.perm_.Height(), 1 );
    invPerm_.Resize( P.invPerm_.Height(), 1 );
    Copy( P.swapDests_, swapDests_.Matrix() );
    Copy( P.swapOrigins_, swapOrigins_.Matrix() );
    Copy( P.perm_, perm_.Matrix() );
    Copy( P.invPerm_, invPerm_.Matrix() );
    
    parity_ = P.parity_;
    staleParity_ = P.staleParity_;
    staleInverse_ = P.staleInverse_;
    staleMeta_ = true;

    return *this;
}

const DistPermutation& DistPermutation::operator=( const DistPermutation& P )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))
    SetGrid( P.grid_ );

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

    colMeta_ = P.colMeta_;
    rowMeta_ = P.rowMeta_;
    staleMeta_ = P.staleMeta_;

    return *this;
}

bool DistPermutation::Parity() const
{
    DEBUG_ONLY(CSE cse("DistPermutation::Parity"))
    if( staleParity_ )
    {
        LogicError("Computation of parity not yet enabled");
        staleParity_ = false;
        return parity_; 
    }
    return parity_;
}

Int DistPermutation::Height() const
{ return size_; }

Int DistPermutation::Width() const
{ return size_; }

bool DistPermutation::IsSwapSequence() const
{ return swapSequence_; }

bool DistPermutation::IsImplicitSwapSequence() const
{ return implicitSwapOrigins_; }

template<typename T>
void DistPermutation::PermuteCols( AbstractDistMatrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::PermuteCols"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_(activeInd,ALL) );
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO: Query/maintain the unordered_map
        const Int align = A.RowAlign();
        mpi::Comm comm = A.RowComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = colMeta_.find( key );
        if( data == colMeta_.end() )
        {
            colMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
            data = colMeta_.find( key );
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, data->second );
    }
}

template<typename T>
void DistPermutation::InversePermuteCols
( AbstractDistMatrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::InversePermuteCols"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_(activeInd,ALL) );
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                ColSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO: Query/maintain the unordered_map
        const Int align = A.RowAlign();
        mpi::Comm comm = A.RowComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = colMeta_.find( key );
        if( data == colMeta_.end() )
        {
            colMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
            data = colMeta_.find( key );
        }
        // TODO: Move El::PermuteCols into this class
        El::PermuteCols( A, data->second, true );
    }
}

template<typename T>
void DistPermutation::PermuteRows( AbstractDistMatrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::PermuteRows"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO: Query/maintain the unordered_map
        const Int align = A.ColAlign();
        mpi::Comm comm = A.ColComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = rowMeta_.find( key );
        if( data == rowMeta_.end() )
        {
            rowMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
            data = rowMeta_.find( key );
        }
        // TODO: Move El::PermuteRows into this class
        El::PermuteRows( A, data->second );
    }
}

template<typename T>
void DistPermutation::InversePermuteRows
( AbstractDistMatrix<T>& A, Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::InversePermuteRows"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO: Query/maintain the unordered_map
        const Int align = A.ColAlign();
        mpi::Comm comm = A.ColComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = rowMeta_.find( key );
        if( data == rowMeta_.end() )
        {
            rowMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
            data = rowMeta_.find( key );
        }
        // TODO: Move El::PermuteRows into this class
        El::PermuteRows( A, data->second, true );
    }
}

template<typename T>
void DistPermutation::PermuteSymmetrically
( UpperOrLower uplo,
  AbstractDistMatrix<T>& A,
  bool conjugate,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::PermuteSymmetrically"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

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
void DistPermutation::InversePermuteSymmetrically
( UpperOrLower uplo,
  AbstractDistMatrix<T>& A,
  bool conjugate,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::InversePermuteSymmetrically"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = origins_STAR_STAR.GetLocal(j,0)+offset;
                const Int dest = dests_STAR_STAR.GetLocal(j,0)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
    }
    else
    {
        if( offset != 0 )
            LogicError
            ("General permutations are not supported with nonzero offsets");

        if( staleMeta_ )
        {
            rowMeta_.clear();
            colMeta_.clear();
            staleMeta_ = false;
        }

        // TODO: Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        LogicError("General symmetric permutations are not yet supported");
    }
}

void DistPermutation::Explicit( AbstractDistMatrix<Int>& P ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::Explicit"))
    P.SetGrid( grid_ );
    // TODO: Faster algorithm 
    Identity( P, size_, size_ );
    PermuteRows( P );
}

#define PROTO(T) \
  template void DistPermutation::PermuteCols \
  ( AbstractDistMatrix<T>& A, \
    Int offset ) const; \
  template void DistPermutation::InversePermuteCols \
  ( AbstractDistMatrix<T>& A, \
    Int offset ) const; \
  template void DistPermutation::PermuteRows \
  ( AbstractDistMatrix<T>& A, \
    Int offset ) const; \
  template void DistPermutation::InversePermuteRows \
  ( AbstractDistMatrix<T>& A, \
    Int offset ) const; \
  template void DistPermutation::PermuteSymmetrically \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<T>& A, \
    bool conjugate, \
    Int offset ) const; \
  template void DistPermutation::InversePermuteSymmetrically \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<T>& A, \
    bool conjugate, \
    Int offset ) const;

#include "El/macros/Instantiate.h"

} // namespace El
