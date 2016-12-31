/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
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
    EL_DEBUG_CSE
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A(i,j) += alpha*ASub(iSub,jSub);
        }
    }
}

// TODO(poulson): Adopt the same approach as GetSubmatrix
template<typename T>
void UpdateSubmatrix
(       AbstractDistMatrix<T>& A,
  const vector<Int>& I,
  const vector<Int>& J,
        T alpha,
  const AbstractDistMatrix<T>& ASub )
{
    EL_DEBUG_CSE
    // TODO(poulson): Intelligently pick the redundant rank to pack from?
    if( ASub.RedundantRank() == 0 )
    {
        const Int ASubLocalHeight = ASub.LocalHeight();
        const Int ASubLocalWidth = ASub.LocalWidth();
        auto& ASubLoc = ASub.LockedMatrix();
        A.Reserve( ASubLocalHeight*ASubLocalWidth );
        for( Int jLoc=0; jLoc<ASubLocalWidth; ++jLoc )
        {
            const Int jSub = ASub.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<ASubLocalHeight; ++iLoc )
            {
                const Int iSub = ASub.GlobalRow(iLoc);
                A.QueueUpdate( I[iSub], J[jSub], alpha*ASubLoc(iLoc,jLoc) );
            }
        }
    }
    A.ProcessQueues();
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void UpdateSubmatrix \
  (       Matrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          T alpha, \
    const Matrix<T>& ASub ); \
  EL_EXTERN template void UpdateSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const vector<Int>& I, \
    const vector<Int>& J, \
          T alpha, \
    const AbstractDistMatrix<T>& ASub );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_UPDATESUBMATRIX_HPP
