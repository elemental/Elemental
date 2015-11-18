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

template<typename T>
void GeneralPurpose
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::GeneralPurpose"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( A.Grid().Size() == 1 )
    {
        B.Resize( height, width );
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        // TODO: Break into smaller pieces to avoid excessive memory usage?
        Zeros( B, height, width );
        if( A.RedundantRank() == 0 )
        {
            B.Reserve( localHeight*localWidth );
            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
                {
                    const Int i = A.GlobalRow(iLoc);
                    B.QueueUpdate( i, j, ABuf[iLoc+jLoc*ALDim] );
                }
            }    
        }

        const bool includeViewers = (A.Grid() != B.Grid());
        B.ProcessQueues( includeViewers );
    }
}

#define PROTO(T) \
  template void GeneralPurpose \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
