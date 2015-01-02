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
void GetSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );

    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            ASub.Set( iSub, jSub, A.Get(i,j) );
        }
    }
}

template<typename T>
void GetRealPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );

    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            ASub.Set( iSub, jSub, A.GetRealPart(i,j) );
        }
    }
}

template<typename T>
void GetImagPartOfSubmatrix
( const Matrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        Matrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );

    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            ASub.Set( iSub, jSub, A.GetImagPart(i,j) );
        }
    }
}

template<typename T>
Matrix<T> GetSubmatrix
( const Matrix<T>& A, const std::vector<Int>& I, const std::vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    Matrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>> GetRealPartOfSubmatrix
( const Matrix<T>& A, const std::vector<Int>& I, const std::vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    Matrix<Base<T>> ASub;
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>> GetImagPartOfSubmatrix
( const Matrix<T>& A, const std::vector<Int>& I, const std::vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    Matrix<Base<T>> ASub;
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<T>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = WriteProxy<T,STAR,STAR>(&ASubPre);
    auto& ASub = *ASubPtr;

    // TODO: Make the following more efficient for non [STAR,STAR]

    ASub.SetGrid( A.Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
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
                        ASub.SetLocal( iSub, jSub, A.GetLocal(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, A.DistComm() );
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, A.Root(), A.CrossComm() );
}

template<typename T>
void GetRealPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = WriteProxy<Base<T>,STAR,STAR>(&ASubPre);
    auto& ASub = *ASubPtr;

    // TODO: Make the following more efficient for non [STAR,STAR]

    ASub.SetGrid( A.Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
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
                        ASub.SetLocal
                        ( iSub, jSub, A.GetLocalRealPart(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, A.DistComm() );
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, A.Root(), A.CrossComm() );
}

template<typename T>
void GetImagPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J, 
        AbstractDistMatrix<Base<T>>& ASubPre )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    auto ASubPtr = WriteProxy<Base<T>,STAR,STAR>(&ASubPre);
    auto& ASub = *ASubPtr;

    // TODO: Make the following more efficient for non [STAR,STAR]

    ASub.SetGrid( A.Grid() );
    ASub.Resize( m, n, m );
    Zeros( ASub, m, n );
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
                        ASub.SetLocal
                        ( iSub, jSub, A.GetLocalImagPart(iLoc,jLoc) );
                    }
                }
            }
        }
        // Sum over the distribution communicator
        mpi::AllReduce( ASub.Buffer(), m*n, A.DistComm() );
    }
    // Broadcast over the cross communicator
    mpi::Broadcast( ASub.Buffer(), m*n, A.Root(), A.CrossComm() );
}

template<typename T>
DistMatrix<T,STAR,STAR> GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J )
{
    DistMatrix<T,STAR,STAR> ASub( A.Grid() );
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J )
{
    DistMatrix<Base<T>,STAR,STAR> ASub( A.Grid() );
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const std::vector<Int>& I, const std::vector<Int>& J )
{
    DistMatrix<Base<T>,STAR,STAR> ASub( A.Grid() );
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

#define PROTO(T) \
  template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          Matrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          Matrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          Matrix<Base<T>>& ASub ); \
  template Matrix<T> GetSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J ); \
  template Matrix<Base<T>> GetRealPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J ); \
  template Matrix<Base<T>> GetImagPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J ); \
  template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          AbstractDistMatrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          AbstractDistMatrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J, \
          AbstractDistMatrix<Base<T>>& ASub ); \
  template DistMatrix<T,STAR,STAR> GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J ); \
  template DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J ); \
  template DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const std::vector<Int>& I, const std::vector<Int>& J );

#include "El/macros/Instantiate.h"

} // namespace El
