/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PERM_DISTPERMUTATION_HPP
#define EL_PERM_DISTPERMUTATION_HPP

#include <map>

namespace El {

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

// TODO(poulson): Convert to accepting Grid rather than mpi::Comm
class DistPermutation
{
public:
    DistPermutation( const Grid& g=Grid::Default() );

    void SetGrid( const Grid& g );

    void Empty();
    void MakeIdentity( Int size );

    void ReserveSwaps( Int maxSwaps );
    void MakeArbitrary() const;

    const DistPermutation& operator=( const Permutation& p );
    const DistPermutation& operator=( const DistPermutation& p );

    void Swap( Int origin, Int dest );
    void SwapSequence( const DistPermutation& perm, Int offset=0 );
    void SwapSequence
    ( const ElementalMatrix<Int>& swapOrigins,
      const ElementalMatrix<Int>& swapDests, Int offset=0 );

    void ImplicitSwapSequence
    ( const ElementalMatrix<Int>& swapDests, Int offset=0 );

    // Explicit (pre)image queries or modifications force the permutation to
    // switch from an explicit sequence of swaps to storing an explicit
    // permutation
    Int Image( Int origin ) const;
    Int Preimage( Int dest ) const;
    void SetImage( Int origin, Int dest );

    // Since local image queries require the permutation to already be explicit,
    // calling this routine when the permutation is a swap sequence will throw
    // an error. These routines are [ADVANCED] since they require knowledge of
    // the distribution of the perm_ and invPerm_ distributed matrices.
    Int LocalImage( Int localOrigin ) const;
    Int LocalPreimage( Int localDest ) const;

    // The following return the same result but follow the usual convention
    Int Height() const;
    Int Width() const;

    bool Parity() const;
    bool IsSwapSequence() const;
    bool IsImplicitSwapSequence() const;

    // NOTE: This is only valid if IsImplicitSwapSequence() is true, otherwise
    //       it is implicitly [0,...,numSwaps-1]
    const DistMatrix<Int,VC,STAR> SwapOrigins() const;
    // NOTE: This is only valid if IsSwapSequence() is true
    const DistMatrix<Int,VC,STAR> SwapDestinations() const;

    template<typename T>
    void PermuteCols
    ( AbstractDistMatrix<T>& A,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteCols
    ( AbstractDistMatrix<T>& A,
      Int offset=0 ) const;

    template<typename T>
    void PermuteRows
    ( AbstractDistMatrix<T>& A,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteRows
    ( AbstractDistMatrix<T>& A,
      Int offset=0 ) const;

    template<typename T>
    void PermuteSymmetrically
    ( UpperOrLower uplo,
      AbstractDistMatrix<T>& A,
      bool conjugate=false,
      Int offset=0 ) const;
    template<typename T>
    void InversePermuteSymmetrically
    ( UpperOrLower uplo,
      AbstractDistMatrix<T>& A,
      bool conjugate=false,
      Int offset=0 ) const;

    // Form the permutation vector p so that P A = A(p,:)
    void ExplicitVector( AbstractDistMatrix<Int>& p ) const;

    // Form the permutation matrix P so that P A = A(p,:)
    void ExplicitMatrix( AbstractDistMatrix<Int>& P ) const;

private:
    Int size_=0;
    const Grid* grid_;
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
    mutable DistMatrix<Int,VC,STAR> swapDests_, swapOrigins_;

    // Only used if swapSequence_=false
    // --------------------------------
    mutable DistMatrix<Int,VC,STAR> perm_;
    mutable DistMatrix<Int,VC,STAR> invPerm_;
    mutable bool staleInverse_=true;

    // Use the alignment and communicator as a key
    typedef std::pair<Int,mpi::Comm> keyType_;
    mutable std::map<keyType_,PermutationMeta> rowMeta_, colMeta_;
    mutable bool staleMeta_=false;
};

} // namespace El

#endif // ifndef EL_PERM_DISTPERMUTATION_HPP
