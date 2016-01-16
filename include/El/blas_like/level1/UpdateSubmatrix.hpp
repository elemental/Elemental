/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_UPDATESUBMATRIX_HPP
#define EL_BLAS_UPDATESUBMATRIX_HPP

namespace El {

template<typename T>
void UpdateSubmatrix
(       Matrix<T>& A, 
  const vector<Int>& I, 
  const vector<Int>& J, 
        T alpha, 
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("UpdateSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.Update( i, j, alpha*ASub.Get(iSub,jSub) );
        }
    }
}

// TODO: Adopt the same approach as GetSubmatrix
template<typename T>
void UpdateSubmatrix
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I,
  const vector<Int>& J, 
        T alpha,
  const AbstractDistMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("UpdateSubmatrix"))
    // TODO: Intelligently pick the redundant rank to pack from?
    if( ASub.RedundantRank() == 0 )
    {
        const Int ASubLocalHeight = ASub.LocalHeight();
        const Int ASubLocalWidth = ASub.LocalWidth();
        A.Reserve( ASubLocalHeight*ASubLocalWidth );
        for( Int jLoc=0; jLoc<ASubLocalWidth; ++jLoc )
        {
            const Int jSub = ASub.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<ASubLocalHeight; ++iLoc )
            {
                const Int iSub = ASub.GlobalRow(iLoc);
                A.QueueUpdate
                ( I[iSub], J[jSub], alpha*ASub.GetLocal(iLoc,jLoc) );
            }
        }
    }
    A.ProcessQueues();
}

} // namespace El

#endif // ifndef EL_BLAS_UPDATESUBMATRIX_HPP
