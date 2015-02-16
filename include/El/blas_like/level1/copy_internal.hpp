/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS1_COPYINTERNALHPP
#define EL_BLAS1_COPYINTERNALHPP

namespace El {

namespace copy {

template<typename T>
void Exchange
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B,
  int sendRank, int recvRank, mpi::Comm comm );

template<typename T,Dist U,Dist V>
void Translate( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );
template<typename T,Dist U,Dist V>
void Translate( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,U,V>& B );

template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );
template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,STAR,STAR>& A, DistMatrix<T,STAR,STAR>& B );
// The fallback case that simply throws an exception
template<typename T,Dist U,Dist V>
void TranslateBetweenGrids
( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );

// NOTE: Only instantiated for (U,V)=(MC,MR) and (U,V)=(MR,MC)
template<typename T,Dist U,Dist V>
void ColwiseVectorExchange
( const DistMatrix<T,ProductDist<U,V>(),STAR>& A,
        DistMatrix<T,ProductDist<V,U>(),STAR>& B );
template<typename T,Dist U,Dist V>
void RowwiseVectorExchange
( const DistMatrix<T,STAR,ProductDist<U,V>()>& A,
        DistMatrix<T,STAR,ProductDist<V,U>()>& B );

// NOTE: Only instantiated for (U,V)=(MC,MR) and (U,V)=(MR,MC)
template<typename T,Dist U,Dist V>
void TransposeDist( const DistMatrix<T,U,V>& A, DistMatrix<T,V,U>& B );

template<typename T,Dist U,Dist V>
void Filter
( const DistMatrix<T,Collect<U>(),Collect<V>()>& A,
        DistMatrix<T,        U,           V   >& B );
template<typename T,Dist U,Dist V>
void Filter
( const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A,
        BlockDistMatrix<T,        U,           V   >& B );

// (V,Collect(U)) |-> (U,V)
template<typename T>
void ColFilter
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B );
template<typename T>
void ColFilter
( const AbstractBlockDistMatrix<T>& A,
        AbstractBlockDistMatrix<T>& B );

// (U,Collect(V)) |-> (U,V)
template<typename T>
void RowFilter
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B );
template<typename T>
void RowFilter
( const AbstractBlockDistMatrix<T>& A,
        AbstractBlockDistMatrix<T>& B );

// (Partial(U),V) |-> (U,V)
template<typename T>
void PartialColFilter
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void PartialColFilter
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

// (U,Partial(V)) |-> (U,V)
template<typename T>
void PartialRowFilter
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void PartialRowFilter
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

template<typename T,Dist U,Dist V>
void AllGather
( const DistMatrix<T,        U,           V   >& A, 
        DistMatrix<T,Collect<U>(),Collect<V>()>& B );
template<typename T,Dist U,Dist V>
void AllGather
( const BlockDistMatrix<T,        U,           V   >& A, 
        BlockDistMatrix<T,Collect<U>(),Collect<V>()>& B );

// (U,V) |-> (Collect(U),V)
template<typename T>
void ColAllGather
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void ColAllGather
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

// (U,V) |-> (U,Collect(V))
template<typename T>
void RowAllGather
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void RowAllGather
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

template<typename T,Dist U,Dist V>
void PartialColAllGather
( const DistMatrix<T,        U,   V>& A,
        DistMatrix<T,Partial<U>(),V>& B );
template<typename T,Dist U,Dist V>
void PartialColAllGather
( const BlockDistMatrix<T,        U,   V>& A,
        BlockDistMatrix<T,Partial<U>(),V>& B );

// (U,V) |-> (U,Partial(V))
template<typename T>
void PartialRowAllGather
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );
template<typename T>
void PartialRowAllGather
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        DistMatrix<T,        U,                     V   >& B );
template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        BlockDistMatrix<T,        U,                     V   >& B );

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A,
        DistMatrix<T,                U,             V   >& B );
template<typename T,Dist U,Dist V>
void RowAllToAllDemote
( const BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A,
        BlockDistMatrix<T,                U,             V   >& B );

template<typename T,Dist U,Dist V>
void ColAllToAllPromote
( const DistMatrix<T,        U,                     V   >& A,
        DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& B );
template<typename T,Dist U,Dist V>
void ColAllToAllPromote
( const BlockDistMatrix<T,        U,                     V   >& A,
        BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& B );

template<typename T,Dist U,Dist V>
void RowAllToAllPromote
( const DistMatrix<T,                U,             V   >& A,
        DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B );
template<typename T,Dist U,Dist V>
void RowAllToAllPromote
( const BlockDistMatrix<T,                U,             V   >& A,
        BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B );

template<typename T>
void Gather
( const AbstractDistMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC>& B );
template<typename T>
void Gather
( const AbstractBlockDistMatrix<T>& A,
        BlockDistMatrix<T,CIRC,CIRC>& B );

template<typename T>
void Scatter
( const DistMatrix<T,CIRC,CIRC>& A,
        AbstractDistMatrix<T>& B );
template<typename T>
void Scatter
( const BlockDistMatrix<T,CIRC,CIRC>& A,
        AbstractBlockDistMatrix<T>& B );

template<typename T>
void Scatter
( const DistMatrix<T,CIRC,CIRC>& A,
        DistMatrix<T,STAR,STAR>& B );
template<typename T>
void Scatter
( const BlockDistMatrix<T,CIRC,CIRC>& A,
        BlockDistMatrix<T,STAR,STAR>& B );

} // namespace copy

} // namespace El

#endif // ifndef EL_BLAS1_COPYINTERNALHPP
