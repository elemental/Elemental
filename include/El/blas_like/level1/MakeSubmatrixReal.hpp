/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_MAKESUBMATRIXREAL_HPP
#define EL_BLAS_MAKESUBMATRIXREAL_HPP

namespace El {

template<typename T>
void MakeSubmatrixReal
(       Matrix<T>& A, 
  const vector<Int>& I,
  const vector<Int>& J )
{
    DEBUG_ONLY(CSE cse("MakeSubmatrixReal"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.MakeReal( i, j );
        }
    }
}

template<typename T>
void MakeSubmatrixReal
(       AbstractDistMatrix<T>& A, 
  const vector<Int>& I,
  const vector<Int>& J )
{
    DEBUG_ONLY(CSE cse("MakeSubmatrixReal"))
    const Int m = I.size();
    const Int n = J.size();

    if( A.Participating() )
    {
        // Fill in our locally-owned entries
        for( Int jSub=0; jSub<n; ++jSub )
        {
            const Int j = J[jSub];
            if( A.IsLocalCol(j) )
            {
                const Int jLoc = A.LocalCol(j);
                for( Int iSub=0; iSub<m; ++iSub )
                {
                    const Int i = I[iSub];
                    if( A.IsLocalRow(i) )
                    {
                        const Int iLoc = A.LocalRow(i);
                        A.MakeLocalReal( iLoc, jLoc );
                    }
                }
            }
        }
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void MakeSubmatrixReal \
  ( Matrix<T>& A, const vector<Int>& I, const vector<Int>& J ); \
  EL_EXTERN template void MakeSubmatrixReal \
  ( AbstractDistMatrix<T>& A, const vector<Int>& I, const vector<Int>& J );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_MAKESUBMATRIXREAL_HPP
