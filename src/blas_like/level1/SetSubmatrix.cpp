/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void SetSubmatrix
(       Matrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("SetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.Set( i, j, ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void SetSubmatrix
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
  const AbstractDistMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("SetSubmatrix"))
    const Grid& g = A.Grid();
    mpi::Comm comm = g.ViewingComm();
    const int commSize = mpi::Size( comm );

    // TODO: Intelligently pick the redundant rank to pack from?

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    if( ASub.RedundantRank() == 0 )
    {
        for( Int jLoc=0; jLoc<ASub.LocalWidth(); ++jLoc )
        {
            const Int jSub = ASub.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<ASub.LocalHeight(); ++iLoc )
            {
                const Int iSub = ASub.GlobalRow(iLoc);
                const int owner = 
                  g.VCToViewing(
                    g.CoordsToVC
                    ( A.ColDist(), A.RowDist(),
                      A.Owner(I[iSub],J[jSub]), A.Root() )
                  );
                ++sendCounts[owner];
            }
        }
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    if( ASub.RedundantRank() == 0 )
    {
        auto offs = sendOffs;
        for( Int jLoc=0; jLoc<ASub.LocalWidth(); ++jLoc )
        {
            const Int jSub = ASub.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<ASub.LocalHeight(); ++iLoc )
            {
                const Int iSub = ASub.GlobalRow(iLoc);
                const int owner = 
                  g.VCToViewing(
                    g.CoordsToVC
                    ( A.ColDist(), A.RowDist(),
                      A.Owner(I[iSub],J[jSub]), A.Root() )
                  );
                const T value = ASub.GetLocal(iLoc,jLoc);
                sendBuf[offs[owner]++] = Entry<T>{ I[iSub], J[jSub], value };
            }
        }
    }

    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll(sendBuf,sendCounts,sendOffs,comm);
    Int recvBufSize = recvBuf.size();
    mpi::Broadcast( recvBufSize, 0, A.RedundantComm() );
    recvBuf.resize( recvBufSize );
    mpi::Broadcast( recvBuf.data(), recvBufSize, 0, A.RedundantComm() );
    for( auto& entry : recvBuf )
        A.Set( entry );
}

// TODO: DistMultiVec version similar to GetSubmatrix implementation

#define PROTO(T) \
  template void SetSubmatrix \
  (       Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
    const Matrix<T>& ASub ); \
  template void SetSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
    const AbstractDistMatrix<T>& ASub );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
