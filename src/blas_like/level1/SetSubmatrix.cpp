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
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("SetSubmatrix"))
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
void SetRealPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("SetRealPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.SetRealPart( i, j, ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void SetImagPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("SetImagPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.SetImagPart( i, j, ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void SetSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<T>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("SetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = ReadProxy<T,STAR,STAR>(&ASubPre);
    const auto& ASub = *ASubPtr;

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
                        A.SetLocal( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void SetRealPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("SetRealPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = ReadProxy<Base<T>,STAR,STAR>(&ASubPre);
    const auto& ASub = *ASubPtr;

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
                        A.SetLocalRealPart
                        ( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void SetImagPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  const AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("SetImagPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = ReadProxy<Base<T>,STAR,STAR>(&ASubPre);
    const auto& ASub = *ASubPtr;

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
                        A.SetLocalImagPart
                        ( iLoc, jLoc, ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

#define PROTO(T) \
  template void SetSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const Matrix<T>& ASub ); \
  template void SetRealPartOfSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const Matrix<Base<T>>& ASub ); \
  template void SetImagPartOfSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const Matrix<Base<T>>& ASub ); \
  template void SetSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const AbstractDistMatrix<T>& ASub ); \
  template void SetRealPartOfSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const AbstractDistMatrix<Base<T>>& ASub ); \
  template void SetImagPartOfSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    const AbstractDistMatrix<Base<T>>& ASub );

#include "El/macros/Instantiate.h"

} // namespace El
