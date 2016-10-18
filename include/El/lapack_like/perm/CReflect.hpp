/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_PERM_CREFLECT_C_HPP
#define EL_LAPACK_PERM_CREFLECT_C_HPP

namespace El {

inline ElPermutationMeta CReflect( const PermutationMeta& meta )
{
    ElPermutationMeta metaC;    

    metaC.align = meta.align;
    metaC.comm = meta.comm.comm;

    const Int commSize = mpi::Size( meta.comm );
    metaC.sendCounts = new int[commSize];
    metaC.sendDispls = new int[commSize];
    metaC.recvCounts = new int[commSize];
    metaC.recvDispls = new int[commSize];
    MemCopy( metaC.sendCounts, meta.sendCounts.data(), commSize );
    MemCopy( metaC.sendDispls, meta.sendDispls.data(), commSize ); 
    MemCopy( metaC.recvCounts, meta.recvCounts.data(), commSize );
    MemCopy( metaC.recvDispls, meta.recvDispls.data(), commSize );

    metaC.numSendIdx = meta.sendIdx.size();
    metaC.numRecvIdx = meta.recvIdx.size();
    metaC.sendIdx   = new int[metaC.numSendIdx];
    metaC.sendRanks = new int[metaC.numSendIdx];
    metaC.recvIdx   = new int[metaC.numRecvIdx];
    metaC.recvRanks = new int[metaC.numRecvIdx];
    MemCopy( metaC.sendIdx,   meta.sendIdx.data(),   metaC.numSendIdx );
    MemCopy( metaC.sendRanks, meta.sendRanks.data(), metaC.numSendIdx );
    MemCopy( metaC.recvIdx,   meta.recvIdx.data(),   metaC.numRecvIdx );
    MemCopy( metaC.recvRanks, meta.recvRanks.data(), metaC.numRecvIdx );

    return metaC;
}

inline PermutationMeta CReflect( const ElPermutationMeta& metaC )
{
    PermutationMeta meta;

    meta.align = metaC.align;
    meta.comm = metaC.comm;

    int commSize;
    MPI_Comm_size( metaC.comm, &commSize );
    meta.sendCounts = 
        vector<int>( metaC.sendCounts, metaC.sendCounts+commSize );
    meta.sendDispls = 
        vector<int>( metaC.sendDispls, metaC.sendDispls+commSize );
    meta.recvCounts =
        vector<int>( metaC.recvCounts, metaC.recvCounts+commSize );
    meta.recvDispls =
        vector<int>( metaC.recvDispls, metaC.recvDispls+commSize );

    meta.sendIdx = 
        vector<int>( metaC.sendIdx, metaC.sendIdx+metaC.numSendIdx );
    meta.sendRanks =
        vector<int>( metaC.sendRanks, metaC.sendRanks+metaC.numSendIdx );
    meta.recvIdx =
        vector<int>( metaC.recvIdx, metaC.recvIdx+metaC.numRecvIdx );
    meta.recvRanks =
        vector<int>( metaC.recvRanks, metaC.recvRanks+metaC.numRecvIdx );

    return meta;
}

inline Permutation* CReflect( ElPermutation p )
{ return EL_RC(Permutation*,p); }
inline ElPermutation CReflect( Permutation* p )
{ return (ElPermutation)EL_RC(struct ElPermutationDummy*,p); }

inline DistPermutation* CReflect( ElDistPermutation p )
{ return EL_RC(DistPermutation*,p); }
inline ElDistPermutation CReflect( DistPermutation* p )
{ return (ElDistPermutation)EL_RC(struct ElDistPermutationDummy*,p); }

inline const Permutation* CReflect( ElConstPermutation p )
{ return EL_RC(const Permutation*,p); }
inline ElConstPermutation CReflect( const Permutation* p )
{ return (ElConstPermutation)EL_RC(const struct ElPermutationDummy*,p); }

inline const DistPermutation* CReflect( ElConstDistPermutation p )
{ return EL_RC(const DistPermutation*,p); }
inline ElConstDistPermutation CReflect( const DistPermutation* p )
{ return (ElConstDistPermutation)
         EL_RC(const struct ElDistPermutationDummy*,p); }

} // namespace El

#endif // ifndef EL_LAPACK_PERM_CREFLECT_C_HPP
