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
void PartialColFilter
( const DistMatrix<T,V,Partial<U>()>& A, 
        DistMatrix<T,U,        V   >& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("transpose::PartialColFilter"))
    DistMatrix<T,V,U> AFilt( A.Grid() );
    if( B.ColConstrained() )
        AFilt.AlignRowsWith( B, false );
    if( B.RowConstrained() )
        AFilt.AlignColsWith( B, false );
    copy::PartialRowFilter( A, AFilt );
    if( !B.ColConstrained() )
        B.AlignColsWith( AFilt, false );
    if( !B.RowConstrained() )
        B.AlignRowsWith( AFilt, false );
    B.Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), B.Matrix(), conjugate );
}

template<typename T,Dist U,Dist V>
void PartialColFilter
( const BlockDistMatrix<T,V,Partial<U>()>& A, 
        BlockDistMatrix<T,U,        V   >& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("transpose::PartialColFilter"))
    BlockDistMatrix<T,V,U> AFilt( A.Grid() );
    if( B.ColConstrained() )
        AFilt.AlignRowsWith( B, false );
    if( B.RowConstrained() )
        AFilt.AlignColsWith( B, false );
    copy::PartialRowFilter( A, AFilt );
    if( !B.ColConstrained() )
        B.AlignColsWith( AFilt, false );
    if( !B.RowConstrained() )
        B.AlignRowsWith( AFilt, false );
    B.Resize( A.Width(), A.Height() );
    Transpose( AFilt.LockedMatrix(), B.Matrix(), conjugate );
}

#define PROTO_DIST(T,U,V) \
  template void PartialColFilter \
  ( const DistMatrix<T,V,Partial<U>()>& A, \
          DistMatrix<T,U,        V   >& B, bool conjugate ); \
  template void PartialColFilter \
  ( const BlockDistMatrix<T,V,Partial<U>()>& A, \
          BlockDistMatrix<T,U,        V   >& B, bool conjugate );

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
