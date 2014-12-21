/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace trans_axpy {

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha, const DistMatrix<T,V,Collect<U>()>& A, 
                 DistMatrix<T,U,        V   >& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("trans_axpy::ColSumScatter"))
    DistMatrix<T,V,U> ASumFilt( A.Grid() );
    if( B.ColConstrained() )
        ASumFilt.AlignRowsWith( B, false );
    if( B.RowConstrained() )
        ASumFilt.AlignColsWith( B, false );
    copy::RowSumScatter( A, ASumFilt );
    if( !B.ColConstrained() )
        B.AlignColsWith( ASumFilt, false );
    if( !B.RowConstrained() )
        B.AlignRowsWith( ASumFilt, false );
    TransposeAxpy( alpha, ASumFilt.LockedMatrix(), B.Matrix(), conjugate );
}

#define PROTO_DIST(T,U,V) \
  template void ColSumScatter \
  ( T alpha, const DistMatrix<T,V,Collect<U>()>& A, \
                   DistMatrix<T,U,        V   >& B, bool conjugate );

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

} // namespace trans_axpy
} // namespace El
