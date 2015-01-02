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
void UpdateSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  T alpha, const Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateSubmatrix"))
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

template<typename T>
void UpdateRealPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateRealPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.UpdateRealPart( i, j, alpha*ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void UpdateImagPartOfSubmatrix
(       Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateImagPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    // Fill in our locally-owned entries
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            A.UpdateImagPart( i, j, alpha*ASub.Get(iSub,jSub) );
        }
    }
}

template<typename T>
void UpdateSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  T alpha, const AbstractDistMatrix<T>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateSubmatrix"))
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
                        A.UpdateLocal
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void UpdateRealPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateRealPartOfSubmatrix"))
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
                        A.UpdateLocalRealPart
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

template<typename T>
void UpdateImagPartOfSubmatrix
(       AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
  Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("UpdateImagPartOfSubmatrix"))
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
                        A.UpdateLocalImagPart
                        ( iLoc, jLoc, alpha*ASub.GetLocal(iSub,jSub) );
                    }
                }
            }
        }
    }
}

#define PROTO(T) \
  template void UpdateSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    T alpha, const Matrix<T>& ASub ); \
  template void UpdateRealPartOfSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    Base<T> alpha, const Matrix<Base<T>>& ASub ); \
  template void UpdateImagPartOfSubmatrix \
  (       Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    Base<T> alpha, const Matrix<Base<T>>& ASub ); \
  template void UpdateSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    T alpha, const AbstractDistMatrix<T>& ASub ); \
  template void UpdateRealPartOfSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub ); \
  template void UpdateImagPartOfSubmatrix \
  (       AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
    Base<T> alpha, const AbstractDistMatrix<Base<T>>& ASub );

#include "El/macros/Instantiate.h"

} // namespace El
