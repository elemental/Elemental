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
(       AbstractDistMatrix<T>& A,
  const PermutationMeta& oldMeta,
  bool inverse=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.RowComm() != oldMeta.comm )
          LogicError("Invalid communicator in metadata");
      if( A.RowAlign() != oldMeta.align )
          LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int localHeight = A.LocalHeight();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localHeight );

    if( inverse )
    {
        // Fill vectors with the send data
        auto offsets = meta.recvDispls;
        const int totalSend = meta.TotalRecv();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.recvIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int jLoc = meta.recvIdx[send];
            const int rank = meta.recvRanks[send];
            MemCopy( &sendData[offsets[rank]], &ABuf[jLoc*ALDim], localHeight );
            offsets[rank] += localHeight;
        }

        // Communicate all pivot rows
        const int totalRecv = meta.TotalSend();
        vector<T> recvData;
        FastResize( recvData, mpi::Pad(totalRecv) );
        mpi::AllToAll
        ( sendData.data(), meta.recvCounts.data(), meta.recvDispls.data(),
          recvData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
          meta.comm );

        // Unpack the recv data
        offsets = meta.sendDispls;
        const int numRecvs = meta.sendIdx.size();
        for( int recv=0; recv<numRecvs; ++recv )
        {
            const int jLoc = meta.sendIdx[recv];
            const int rank = meta.sendRanks[recv];
            MemCopy( &ABuf[jLoc*ALDim], &recvData[offsets[rank]], localHeight );
            offsets[rank] += localHeight;
        }
    }
    else
    {
        // Fill vectors with the send data
        auto offsets = meta.sendDispls;
        const int totalSend = meta.TotalSend();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.sendIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int jLoc = meta.sendIdx[send];
            const int rank = meta.sendRanks[send];
            MemCopy( &sendData[offsets[rank]], &ABuf[jLoc*ALDim], localHeight );
            offsets[rank] += localHeight;
        }

        // Communicate all pivot rows
        const int totalRecv = meta.TotalRecv();
        vector<T> recvData;
        FastResize( recvData, mpi::Pad(totalRecv) );
        mpi::AllToAll
        ( sendData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
          recvData.data(), meta.recvCounts.data(), meta.recvDispls.data(),
          meta.comm );

        // Unpack the recv data
        offsets = meta.recvDispls;
        const int numRecvs = meta.recvIdx.size();
        for( int recv=0; recv<numRecvs; ++recv )
        {
            const int jLoc = meta.recvIdx[recv];
            const int rank = meta.recvRanks[recv];
            MemCopy( &ABuf[jLoc*ALDim], &recvData[offsets[rank]], localHeight );
            offsets[rank] += localHeight;
        }
    }
}

template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const PermutationMeta& oldMeta,
  bool inverse=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.ColComm() != oldMeta.comm )
          LogicError("Invalid communicator in metadata");
      if( A.ColAlign() != oldMeta.align )
          LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int localWidth = A.LocalWidth();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localWidth );

    if( inverse )
    {
        // Fill vectors with the send data
        auto offsets = meta.recvDispls;
        const int totalSend = meta.TotalRecv();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.recvIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int iLoc = meta.recvIdx[send];
            const int rank = meta.recvRanks[send];

            StridedMemCopy
            ( &sendData[offsets[rank]], 1, &ABuf[iLoc], ALDim, localWidth );
            offsets[rank] += localWidth;
        }

        // Communicate all pivot rows
        const int totalRecv = meta.TotalSend();
        vector<T> recvData;
        FastResize( recvData, mpi::Pad(totalRecv) );
        mpi::AllToAll
        ( sendData.data(), meta.recvCounts.data(), meta.recvDispls.data(),
          recvData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
          meta.comm );

        // Unpack the recv data
        offsets = meta.sendDispls;
        const int numRecvs = meta.sendIdx.size();
        for( int recv=0; recv<numRecvs; ++recv )
        {
            const int iLoc = meta.sendIdx[recv];
            const int rank = meta.sendRanks[recv];
            StridedMemCopy
            ( &ABuf[iLoc], ALDim, &recvData[offsets[rank]], 1,localWidth );
            offsets[rank] += localWidth;
        }
    }
    else
    {
        // Fill vectors with the send data
        auto offsets = meta.sendDispls;
        const int totalSend = meta.TotalSend();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.sendIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int iLoc = meta.sendIdx[send];
            const int rank = meta.sendRanks[send];

            StridedMemCopy
            ( &sendData[offsets[rank]], 1, &ABuf[iLoc], ALDim, localWidth );
            offsets[rank] += localWidth;
        }

        // Communicate all pivot rows
        const int totalRecv = meta.TotalRecv();
        vector<T> recvData;
        FastResize( recvData, mpi::Pad(totalRecv) );
        mpi::AllToAll
        ( sendData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
          recvData.data(), meta.recvCounts.data(), meta.recvDispls.data(),
          meta.comm );

        // Unpack the recv data
        offsets = meta.recvDispls;
        const int numRecvs = meta.recvIdx.size();
        for( int recv=0; recv<numRecvs; ++recv )
        {
            const int iLoc = meta.recvIdx[recv];
            const int rank = meta.recvRanks[recv];
            StridedMemCopy
            ( &ABuf[iLoc], ALDim, &recvData[offsets[rank]], 1,localWidth );
            offsets[rank] += localWidth;
        }
    }
}

void InvertPermutation
( const AbstractDistMatrix<Int>& pPre,
        AbstractDistMatrix<Int>& pInvPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( pPre.Width() != 1 )
          LogicError("p must be a column vector");
    )

    const Int n = pPre.Height();
    pInvPre.Resize( n, 1 );
    if( n == 0 )
        return;

    DistMatrixReadProxy<Int,Int,VC,STAR> pProx( pPre );
    DistMatrixWriteProxy<Int,Int,VC,STAR> pInvProx( pInvPre );
    auto& p = pProx.GetLocked();
    auto& pInv = pInvProx.Get();
    auto& pLoc = p.LockedMatrix();
    auto& pInvLoc = pInv.Matrix();

    EL_DEBUG_ONLY(
      // This is obviously necessary but not sufficient for 'p' to contain
      // a reordering of (0,1,...,n-1).
      const Int range = MaxNorm( p ) + 1;
      if( range != n )
          LogicError("Invalid permutation range");
    )

    // Compute the send counts
    const mpi::Comm colComm = p.ColComm();
    const Int commSize = mpi::Size( colComm );
    vector<int> sendSizes(commSize,0), recvSizes(commSize,0);
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int iDest = pLoc(iLoc);
        const int owner = pInv.RowOwner(iDest);
        sendSizes[owner] += 2; // we'll send the global index and the value
    }
    // Perform a small AllToAll to get the receive counts
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, colComm );
    vector<int> sendOffs, recvOffs;
    const int sendTotal = Scan( sendSizes, sendOffs );
    const int recvTotal = Scan( recvSizes, recvOffs );

    // Pack the send data
    vector<Int> sendBuf(sendTotal);
    auto offsets = sendOffs;
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int i     = p.GlobalRow(iLoc);
        const Int iDest = pLoc(iLoc);
        const int owner = pInv.RowOwner(iDest);
        sendBuf[offsets[owner]++] = iDest;
        sendBuf[offsets[owner]++] = i;
    }

    // Perform the actual exchange
    vector<Int> recvBuf(recvTotal);
    mpi::AllToAll
    ( sendBuf.data(), sendSizes.data(), sendOffs.data(),
      recvBuf.data(), recvSizes.data(), recvOffs.data(), colComm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the received data
    for( Int k=0; k<recvTotal/2; ++k )
    {
        const Int iDest = recvBuf[2*k+0];
        const Int i     = recvBuf[2*k+1];

        const Int iDestLoc = pInv.LocalRow(iDest);
        pInvLoc(iDestLoc) = i;
    }
}

} // anonymous namespace

DistPermutation::DistPermutation( const Grid& g )
: grid_(&g)
{
    EL_DEBUG_CSE
    swapDests_.SetGrid( g );
    swapOrigins_.SetGrid( g );
    perm_.SetGrid( g );
    invPerm_.SetGrid( g );
}

void DistPermutation::SetGrid( const Grid& g )
{
    EL_DEBUG_CSE
    Empty();
    grid_ = &g;
    swapDests_.SetGrid( g );
    swapOrigins_.SetGrid( g );
    perm_.SetGrid( g );
    invPerm_.SetGrid( g );
}

void DistPermutation::Empty()
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

    rowMeta_.clear();
    colMeta_.clear();
    staleMeta_ = false;
}

void DistPermutation::MakeIdentity( Int size )
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

void DistPermutation::ReserveSwaps( Int maxSwaps )
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

void DistPermutation::Swap( Int origin, Int dest )
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

        swapOrigins_.SetGrid( *grid_ );
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

void DistPermutation::SwapSequence( const DistPermutation& P, Int offset )
{
    EL_DEBUG_CSE
    if( P.swapSequence_ )
    {
        const Int numSwapAppends = P.numSwaps_;
        auto activeInd = IR(0,numSwapAppends);

        DistMatrix<Int,STAR,STAR> swapDests_STAR_STAR =
          P.swapDests_(activeInd,ALL);
        auto& swapDestsLoc = swapDests_STAR_STAR.Matrix();

        if( P.implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwapAppends; ++j )
                Swap( j+offset, swapDestsLoc(j)+offset );
        }
        else
        {
            DistMatrix<Int,STAR,STAR> swapOrigins_STAR_STAR =
              P.swapOrigins_(activeInd,ALL);
            auto& swapOriginsLoc = swapOrigins_STAR_STAR.Matrix();
            for( Int j=0; j<numSwapAppends; ++j )
                Swap( swapOriginsLoc(j)+offset, swapDestsLoc(j)+offset );
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

void DistPermutation::SwapSequence
( const ElementalMatrix<Int>& swapOriginsPre,
  const ElementalMatrix<Int>& swapDestsPre,
  Int offset )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Int,Int,STAR,STAR>
      swapOriginsProx( swapOriginsPre ),
      swapDestsProx( swapDestsPre );
    auto& swapOrigins = swapOriginsProx.GetLocked();
    auto& swapDests = swapDestsProx.GetLocked();
    auto& swapOriginsLoc = swapOrigins.LockedMatrix();
    auto& swapDestsLoc = swapDests.LockedMatrix();

    // TODO(poulson): Assert swapOrigins and swapDests are column vectors of
    // same size
    const Int numSwaps = swapDests.Height();
    for( Int k=0; k<numSwaps; ++k )
        Swap( swapOriginsLoc(k)+offset, swapDestsLoc(k)+offset );
}

void DistPermutation::ImplicitSwapSequence
( const ElementalMatrix<Int>& swapDestsPre,
  Int offset )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Int,Int,STAR,STAR> swapDestsProx( swapDestsPre );
    auto& swapDests = swapDestsProx.GetLocked();
    auto& swapDestsLoc = swapDests.LockedMatrix();

    const Int numPrevSwaps = numSwaps_;

    // TODO(poulson): Assert swapOrigins and swapDests are column vectors of
    // same size
    const Int numSwaps = swapDests.Height();
    for( Int k=0; k<numSwaps; ++k )
        Swap( numPrevSwaps+k, swapDestsLoc(k)+offset );
}

Int DistPermutation::LocalImage( Int localOrigin ) const
{
    EL_DEBUG_CSE
    if( swapSequence_ )
        LogicError("Cannot query the local image of a swap sequence");
    if( staleInverse_ )
        LogicError("Cannot query the local image when the inverse is stale");
    return invPerm_.GetLocal( localOrigin, 0 );
}

Int DistPermutation::LocalPreimage( Int localDest ) const
{
    EL_DEBUG_CSE
    if( swapSequence_ )
        LogicError("Cannot query the local preimage of a swap sequence");
    return perm_.GetLocal( localDest, 0 );
}

Int DistPermutation::Image( Int origin ) const
{
    EL_DEBUG_CSE
    MakeArbitrary();
    if( staleInverse_ )
    {
        El::InvertPermutation( perm_, invPerm_ );
        staleInverse_ = false;
    }
    return invPerm_.Get( origin, 0 );
}

Int DistPermutation::Preimage( Int dest ) const
{
    EL_DEBUG_CSE
    MakeArbitrary();
    return perm_.Get( dest, 0 );
}

void DistPermutation::SetImage( Int origin, Int dest )
{
    EL_DEBUG_CSE
    MakeArbitrary();
    perm_.Set( dest, 0, origin );
    invPerm_.Set( origin, 0, dest );
    staleParity_ = true;
    staleMeta_ = true;
}

void DistPermutation::MakeArbitrary() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ )
        return;

    // Compose the swaps into an explicit permutation vector
    // -----------------------------------------------------
    perm_.Resize( size_, 1 );
    auto& permLoc = perm_.Matrix();
    const Int localHeight = perm_.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        permLoc(iLoc) = perm_.GlobalRow(iLoc);
    PermuteRows( perm_ );

    invPerm_.Resize( size_, 1 );
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
    EL_DEBUG_CSE
    if( grid_->Size() != 1 )
        LogicError
        ("Invalid grid size of ",grid_->Size()," for sequential copy");

    size_ = P.size_;

    swapSequence_ = P.swapSequence_;
    numSwaps_ = P.numSwaps_;
    implicitSwapOrigins_ = P.implicitSwapOrigins_;

    swapDests_.Resize( P.swapDests_.Height(), P.swapDests_.Width() );
    Copy( P.swapDests_, swapDests_.Matrix() );

    swapOrigins_.Resize( P.swapOrigins_.Height(), P.swapOrigins_.Width() );
    Copy( P.swapOrigins_, swapOrigins_.Matrix() );

    perm_.Resize( P.perm_.Height(), P.perm_.Width() );
    Copy( P.perm_, perm_.Matrix() );

    invPerm_.Resize( P.invPerm_.Height(), P.invPerm_.Width() );
    Copy( P.invPerm_, invPerm_.Matrix() );

    parity_ = P.parity_;
    staleParity_ = P.staleParity_;
    staleInverse_ = P.staleInverse_;
    staleMeta_ = true;

    return *this;
}

const DistPermutation& DistPermutation::operator=( const DistPermutation& P )
{
    EL_DEBUG_CSE
    SetGrid( *P.grid_ );

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

    colMeta_ = P.colMeta_;
    rowMeta_ = P.rowMeta_;
    staleMeta_ = P.staleMeta_;

    return *this;
}

bool DistPermutation::Parity() const
{
    EL_DEBUG_CSE
    if( staleParity_ )
    {
        if( swapSequence_ )
            LogicError("Unexpected stale parity for a swap sequence");

        // Walking through the process of LU factorization with partial
        // pivoting for a permutation matrix, which never requires a
        // Schur-complement update, yields an algorithm for expressing the
        // inverse of a permutation in terms of a sequence of transpositions in
        // linear time. Note that performing the swaps requires access to the
        // inverse permutation, which can be formed in linear time.
        if( staleInverse_ )
        {
            El::InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        DistMatrix<Int,STAR,STAR> permCopy(perm_), invPermCopy(invPerm_);
        auto& permCopyLoc = permCopy.Matrix();
        auto& invPermCopyLoc = invPermCopy.Matrix();

        parity_ = false;
        for( Int k=0; k<size_; ++k )
        {
            const Int permVal = permCopyLoc(k);
            if( permVal != k )
            {
                parity_ = !parity_;
                const Int invPermVal = invPermCopyLoc(k);
                // We only need to perform half of the swaps
                //      perm[k] <-> perm[invPerm[k]]
                //   invPerm[k] <-> invPerm[perk[k]]
                // since we will not need to access perm[k] and invPerm[k] again
                permCopyLoc(invPermVal) = permVal;
                invPermCopyLoc(permVal) = invPermVal;
            }
        }

        staleParity_ = false;
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

const DistMatrix<Int,VC,STAR> DistPermutation::SwapOrigins() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ || !implicitSwapOrigins_ )
        LogicError("Swap origins are not explicitly stored");
    return swapOrigins_(IR(0,numSwaps_),ALL);
}

const DistMatrix<Int,VC,STAR> DistPermutation::SwapDestinations() const
{
    EL_DEBUG_CSE
    if( !swapSequence_ )
        LogicError("Swap destinations are not explicitly stored");
    return swapDests_(IR(0,numSwaps_),ALL);
}

template<typename T>
void DistPermutation::PermuteCols( AbstractDistMatrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        // TODO(poulson): Introduce an std::map for caching the pivots this
        // process needs to care about to avoid redundant [STAR,STAR] formations
        //
        // TODO(poulson): Decide on the data structure; an std::vector of pivot
        // origins and destinations?

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_(activeInd,ALL) );
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO(poulson): Query/maintain the unordered_map
        const Int align = A.RowAlign();
        mpi::Comm comm = A.RowComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = colMeta_.find( key );
        if( data == colMeta_.end() )
        {
// TODO(poulson): Enable this branch; it apparently is not possible with
// GCC 4.7.1
#ifdef EL_HAVE_STD_EMPLACE
            colMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
#else
            auto newPair =
              std::make_pair(key,PermutationMeta(perm_,invPerm_,align,comm));
            colMeta_.insert( newPair );
#endif
            data = colMeta_.find( key );
        }
        // TODO(poulson): Move El::PermuteCols into this class
        El::PermuteCols( A, data->second );
    }
}

template<typename T>
void DistPermutation::InversePermuteCols
( AbstractDistMatrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR( swapDests_(activeInd,ALL) );
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                ColSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO(poulson): Query/maintain the unordered_map
        const Int align = A.RowAlign();
        mpi::Comm comm = A.RowComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = colMeta_.find( key );
        if( data == colMeta_.end() )
        {
// TODO(poulson): Enable this branch; it apparently is not possible with
// GCC 4.7.1
#ifdef EL_HAVE_STD_EMPLACE
            colMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
#else
            auto newPair =
              std::make_pair(key,PermutationMeta(perm_,invPerm_,align,comm));
            colMeta_.insert( newPair );
#endif
            data = colMeta_.find( key );
        }
        // TODO(poulson): Move El::PermuteCols into this class
        El::PermuteCols( A, data->second, true );
    }
}

template<typename T>
void DistPermutation::PermuteRows( AbstractDistMatrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO(poulson): Query/maintain the unordered_map
        const Int align = A.ColAlign();
        mpi::Comm comm = A.ColComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = rowMeta_.find( key );
        if( data == rowMeta_.end() )
        {
// TODO(poulson): Enable this branch; it apparently is not possible with
// GCC 4.7.1
#ifdef EL_HAVE_STD_EMPLACE
            rowMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
#else
            auto newPair =
              std::make_pair(key,PermutationMeta(perm_,invPerm_,align,comm));
            rowMeta_.insert( newPair );
#endif
            data = rowMeta_.find( key );
        }
        // TODO(poulson): Move El::PermuteRows into this class
        El::PermuteRows( A, data->second );
    }
}

template<typename T>
void DistPermutation::InversePermuteRows
( AbstractDistMatrix<T>& A, Int offset ) const
{
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                El::RowSwap( A, origin, dest );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }

        // TODO(poulson): Query/maintain the unordered_map
        const Int align = A.ColAlign();
        mpi::Comm comm = A.ColComm();
        keyType_ key = std::pair<Int,mpi::Comm>(align,comm);
        auto data = rowMeta_.find( key );
        if( data == rowMeta_.end() )
        {
// TODO(poulson): Enable this branch; it apparently is not possible with
// GCC 4.7.1
#ifdef EL_HAVE_STD_EMPLACE
            rowMeta_.emplace
            ( std::piecewise_construct,
              std::forward_as_tuple(key),
              std::forward_as_tuple(perm_,invPerm_,align,comm) );
#else
            auto newPair =
              std::make_pair(key,PermutationMeta(perm_,invPerm_,align,comm));
            rowMeta_.insert( newPair );
#endif
            data = rowMeta_.find( key );
        }
        // TODO(poulson): Move El::PermuteRows into this class
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
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=0; j<numSwaps_; ++j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
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
    EL_DEBUG_CSE
    // TODO(poulson): Use an (MC,MR) proxy for A?
    if( swapSequence_ )
    {
        const Int height = A.Height();
        const Int width = A.Width();
        if( height == 0 || width == 0 )
            return;

        auto activeInd = IR(0,numSwaps_);

        DistMatrix<Int,STAR,STAR> dests_STAR_STAR = swapDests_(activeInd,ALL);
        auto& destsLoc = dests_STAR_STAR.Matrix();
        if( implicitSwapOrigins_ )
        {
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = j+offset;
                const Int dest = destsLoc(j)+offset;
                SymmetricSwap( uplo, A, origin, dest, conjugate );
            }
        }
        else
        {
            DistMatrix<Int,STAR,STAR> origins_STAR_STAR =
              swapOrigins_(activeInd,ALL);
            auto& originsLoc = origins_STAR_STAR.Matrix();
            for( Int j=numSwaps_-1; j>=0; --j )
            {
                const Int origin = originsLoc(j)+offset;
                const Int dest = destsLoc(j)+offset;
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

        // TODO(poulson): Move El::InversePermutation into this class
        if( staleInverse_ )
        {
            InvertPermutation( perm_, invPerm_ );
            staleInverse_ = false;
        }
        LogicError("General symmetric permutations are not yet supported");
    }
}

void DistPermutation::ExplicitVector( AbstractDistMatrix<Int>& p ) const
{
    EL_DEBUG_CSE
    p.SetGrid( *grid_ );
    if( swapSequence_ )
    {
        p.Resize( size_, 1 );
        auto& pLoc = p.Matrix();
        const Int localHeight = p.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            pLoc(iLoc) = p.GlobalRow(iLoc);
        PermuteRows( p );
    }
    else
    {
        Copy( perm_, p );
    }
}

void DistPermutation::ExplicitMatrix( AbstractDistMatrix<Int>& P ) const
{
    EL_DEBUG_CSE
    P.SetGrid( *grid_ );

    DistMatrix<Int,VC,STAR> p_VC_STAR(*grid_);
    ExplicitVector( p_VC_STAR );
    DistMatrix<Int,STAR,STAR> p( p_VC_STAR );
    auto& pLoc = p.LockedMatrix();

    Zeros( P, size_, size_ );
    for( Int i=0; i<size_; ++i )
        P.Set( i, pLoc(i), 1 );
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

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
