/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El/blas_like/level1/copy_internal.hpp"

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void ColwiseVectorExchange
( const DistMatrix<T,ProductDist<U,V>(),STAR>& A,
        DistMatrix<T,ProductDist<V,U>(),STAR>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColwiseVectorExchange"))
    AssertSameGrids( A, B );
    if( !B.Participating() )
        return;

    const Int distSize = A.DistSize();
    const Int colDiff = A.ColShift() - B.ColShift();
    const Int sendRankB = Mod( B.DistRank()+colDiff, distSize );
    const Int recvRankA = Mod( A.DistRank()-colDiff, distSize );
    const Int recvRankB =
      (recvRankA/A.PartialColStride())+
      (recvRankA%A.PartialColStride())*A.PartialUnionColStride();
    copy::Exchange( A, B, sendRankB, recvRankB, B.DistComm() );
}

#define PROTO_DIST(T,U,V) \
  template void ColwiseVectorExchange<T,U,V> \
  ( const DistMatrix<T,ProductDist<U,V>(),STAR>& A, \
          DistMatrix<T,ProductDist<V,U>(),STAR>& B );

#define PROTO(T) \
  PROTO_DIST(T,MC,MR) \
  PROTO_DIST(T,MR,MC)

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
