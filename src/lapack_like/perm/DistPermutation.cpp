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

    rowMeta_.clear();
    colMeta_.clear();
    staleMeta_ = false;
}

void DistPermutation::ReserveSwaps( Int numSwaps )
{
    DEBUG_ONLY(CSE cse("DistPermutation::ReserveSwaps"))
    Empty();

    swapDests_.SetGrid( grid_ );
    swapDests_.Resize( numSwaps, 1 );
    const Int localHeight = swapDests_.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        swapDests_.SetLocal( iLoc, 0, swapDests_.GlobalRow(iLoc) );

    staleMeta_ = true;
}

void DistPermutation::AppendSwap( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("DistPermutation::AppendSwap"))
    if( !swapSequence_ )
        LogicError("Appending swaps to general permutations not yet enabled");
    if( !staleParity_ && origin != dest )
        parity_ = !parity_;
    staleMeta_ = true;

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

        swapOrigins_.SetGrid( grid_ );
        swapOrigins_.Resize( numSwaps, 1 );
        const Int localHeight = swapOrigins_.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            swapOrigins_.SetLocal( iLoc, 0, swapOrigins_.GlobalRow(iLoc) );

        swapOrigins_.Set( nextSwapIndex_, 0, origin );
        swapDests_.Set( nextSwapIndex_, 0, dest );
        ++nextSwapIndex_;
    }
}

void DistPermutation::AppendSwapSequence
( const DistPermutation& perm, Int offset )
{
    DEBUG_ONLY(CSE cse("DistPermutation::AppendSwapSequence"))
    if( !swapSequence_ || !perm.swapSequence_ )
        LogicError("Appending swaps to general permutations not yet enabled");
    
    const Int numSwapAppends = perm.swapDests_.Height();
    DistMatrix<Int,STAR,STAR> swapDests_STAR_STAR = perm.swapDests_;
    if( perm.implicitSwapOrigins_ )
    {
        for( Int j=0; j<numSwapAppends; ++j )
            AppendSwap
            ( nextSwapIndex_+j, swapDests_STAR_STAR.GetLocal(j,0)-offset );
    }
    else
    {
        DistMatrix<Int,STAR,STAR> swapOrigins_STAR_STAR = perm.swapOrigins_;
        for( Int j=0; j<numSwapAppends; ++j )
            AppendSwap
            ( swapOrigins_STAR_STAR.GetLocal(j,0)-offset,
              swapDests_STAR_STAR.GetLocal(j,0)-offset );
    }
    nextSwapIndex_ += numSwapAppends;
}

void DistPermutation::SetImage( Int origin, Int dest )
{
    DEBUG_ONLY(CSE cse("DistPermutation::SetImage"))
    if( swapSequence_ )
        LogicError("Cannot explicitly set the image of a swap sequence");
    perm_.Set( origin, 0, dest );
    invPerm_.Set( dest, 0, origin );
    staleParity_ = true;
    staleMeta_ = true;
}

void DistPermutation::MakeArbitrary( Int domainSize )
{
    DEBUG_ONLY(CSE cse("DistPermutation::MakeArbitrary"))
    if( swapSequence_ )
    {
        Empty();
        swapSequence_ = false;
    }
    perm_.SetGrid( grid_ );
    invPerm_.SetGrid( grid_ );
    perm_.Resize( domainSize, 1 );
    invPerm_.Resize( domainSize, 1 );
    const Int localHeight = perm_.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = perm_.GlobalRow(iLoc);
        perm_.SetLocal( iLoc, 0, i );
        invPerm_.SetLocal( iLoc, 0, i );
    }
    parity_ = false;
    staleParity_ = false;
    staleInverse_ = false;
    staleMeta_ = true;
}

const DistPermutation& DistPermutation::operator=( const Permutation& p )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))
    if( grid_.Size() != 1 )
        LogicError("Invalid grid size for sequential copy");

    swapSequence_ = p.swapSequence_;
    nextSwapIndex_ = p.nextSwapIndex_;
    implicitSwapOrigins_ = p.implicitSwapOrigins_;

    swapDests_.Resize( p.swapDests_.Height(), 1 );
    swapOrigins_.Resize( p.swapOrigins_.Height(), 1 );
    perm_.Resize( p.perm_.Height(), 1 );
    invPerm_.Resize( p.invPerm_.Height(), 1 );
    Copy( p.swapDests_, swapDests_.Matrix() );
    Copy( p.swapOrigins_, swapOrigins_.Matrix() );
    Copy( p.perm_, perm_.Matrix() );
    Copy( p.invPerm_, invPerm_.Matrix() );
    
    parity_ = p.parity_;
    staleParity_ = p.staleParity_;
    staleInverse_ = p.staleInverse_;
    staleMeta_ = true;

    return *this;
}

const DistPermutation& DistPermutation::operator=( const DistPermutation& p )
{
    DEBUG_ONLY(CSE cse("DistPermutation::operator="))
    SetGrid( p.grid_ );

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

    colMeta_ = p.colMeta_;
    rowMeta_ = p.rowMeta_;
    staleMeta_ = p.staleMeta_;

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

bool DistPermutation::IsSwapSequence() const
{ return swapSequence_; }

bool DistPermutation::IsImplicitSwapSequence() const
{ return implicitSwapOrigins_; }

template<typename T>
void DistPermutation::PermuteCols
( AbstractDistMatrix<T>& A,
  bool inverse,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::PermuteCols"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        const Int numSwaps = swapDests_.Height();
        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_ );
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
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
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    ColSwap( A, origin, dest );
                }
            }
        }
    }
    else
    {
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
        El::PermuteCols( A, data->second, inverse );
    }
}

template<typename T>
void DistPermutation::PermuteRows
( AbstractDistMatrix<T>& A,
  bool inverse,
  Int offset ) const
{
    DEBUG_ONLY(CSE cse("DistPermutation::PermuteRows"))
    // TODO: Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        const Int numSwaps = swapDests_.Height();
        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_ );
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
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
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    RowSwap( A, origin, dest );
                }
            }
        }
    }
    else
    {
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
        El::PermuteRows( A, data->second, inverse );
    }
}

template<typename T>
void DistPermutation::PermuteSymmetrically
( UpperOrLower uplo,
  AbstractDistMatrix<T>& A,
  bool conjugate,
  bool inverse,
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

        const Int numSwaps = swapDests_.Height();
        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_ );
        if( inverse )
        {
            if( implicitSwapOrigins_ )
            {
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = j;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=numSwaps-1; j>=0; --j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
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
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }
            else
            {
                DistMatrix<Int,STAR,STAR> origins_STAR_STAR( swapOrigins_ );
                for( Int j=0; j<numSwaps; ++j )
                {
                    const Int origin = origins_STAR_STAR.GetLocal(j,0)-offset;
                    const Int dest = dests_STAR_STAR.GetLocal(j,0)-offset;
                    SymmetricSwap( uplo, A, origin, dest, conjugate );
                }
            }
        }
    }
    else
    {
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
    if( swapSequence_ )
    {
        DistMatrix<Int,VC,STAR> perm(grid_);
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
  template void DistPermutation::PermuteCols \
  ( AbstractDistMatrix<T>& A, \
    bool inverse, \
    Int offset ) const; \
  template void DistPermutation::PermuteRows \
  ( AbstractDistMatrix<T>& A, \
    bool inverse, \
    Int offset ) const; \
  template void DistPermutation::PermuteSymmetrically \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<T>& A, \
    bool conjugate, \
    bool inverse, \
    Int offset ) const;

#include "El/macros/Instantiate.h"

} // namespace El
