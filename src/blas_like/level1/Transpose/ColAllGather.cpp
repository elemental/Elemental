/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace transpose {

template<typename T,Dist U,Dist V>
void ColAllGather
( const DistMatrix<T,U,V           >& A, 
        DistMatrix<T,V,Collect<U>()>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("transpose::ColAllGather"))
    DistMatrix<T,V,U> ATrans( A.Grid() );
    ATrans.AlignWith( A );
    ATrans.Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans.Matrix(), conjugate );
    copy::RowAllGather( ATrans, B );
}

template<typename T,Dist U,Dist V>
void ColAllGather
( const BlockDistMatrix<T,U,V           >& A, 
        BlockDistMatrix<T,V,Collect<U>()>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("transpose::ColAllGather"))
    BlockDistMatrix<T,V,U> ATrans( A.Grid() );
    ATrans.AlignWith( A );
    ATrans.Resize( A.Width(), A.Height() );
    Transpose( A.LockedMatrix(), ATrans.Matrix(), conjugate );
    copy::RowAllGather( ATrans, B );
}

#define PROTO_DIST(T,U,V) \
  template void ColAllGather \
  ( const DistMatrix<T,U,V           >& A, \
          DistMatrix<T,V,Collect<U>()>& B, bool conjugate ); \
  template void ColAllGather \
  ( const BlockDistMatrix<T,U,V           >& A, \
          BlockDistMatrix<T,V,Collect<U>()>& B, bool conjugate );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR)

#include "El/macros/Instantiate.h"

} // namespace transpose
} // namespace El
