/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_CONCATENATE_HPP
#define EL_BLAS_CONCATENATE_HPP

namespace El {

template<typename T>
void HCat
( const Matrix<T>& A,
  const Matrix<T>& B,
        Matrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");
    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    C.Resize( m, nA+nB );
    Zero( C );
    auto CL = C( IR(0,m), IR(0,nA)     );
    auto CR = C( IR(0,m), IR(nA,nA+nB) );
    CL = A;
    CR = B;
}

template<typename T>
void VCat
( const Matrix<T>& A,
  const Matrix<T>& B,
        Matrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    C.Resize( mA+mB, n );
    Zero( C );
    auto CT = C( IR(0,mA),     IR(0,n) );
    auto CB = C( IR(mA,mA+mB), IR(0,n) );
    CT = A;
    CB = B;
}

template<typename T>
inline void HCat
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& CPre )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");
    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    DistMatrixWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    C.Resize( m, nA+nB );
    Zero( C );
    auto CL = C( IR(0,m), IR(0,nA)     );
    auto CR = C( IR(0,m), IR(nA,nA+nB) );
    CL = A;
    CR = B;
}

template<typename T>
void VCat
( const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& CPre )
{
    EL_DEBUG_CSE
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    DistMatrixWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    C.Resize( mA+mB, n );
    Zero( C );
    auto CT = C( IR(0,mA),     IR(0,n) );
    auto CB = C( IR(mA,mA+mB), IR(0,n) );
    CT = A;
    CB = B;
}

template<typename T>
void HCat
( const SparseMatrix<T>& A,
  const SparseMatrix<T>& B,
        SparseMatrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");

    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    C.Resize( m, nA+nB );
    Zero( C );
    C.Reserve( numEntriesA+numEntriesB );
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueUpdate( B.Row(e), B.Col(e)+nA, B.Value(e) );
    C.ProcessQueues();
}

template<typename T>
void VCat
( const SparseMatrix<T>& A,
  const SparseMatrix<T>& B,
        SparseMatrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    C.Resize( mA+mB, n );
    Zero( C );
    C.Reserve( numEntriesA+numEntriesB );
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueUpdate( B.Row(e)+mA, B.Col(e), B.Value(e) );
    C.ProcessQueues();
}

template<typename T>
void HCat
( const DistSparseMatrix<T>& A,
  const DistSparseMatrix<T>& B,
        DistSparseMatrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");
    /*
    if( A.Grid().Comm() != B.Grid().Comm() )
        LogicError("A and B had different communicators");
    */

    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesB = B.NumLocalEntries();
    C.SetGrid( A.Grid() );
    C.Resize( m, nA+nB );
    Zero( C );
    C.Reserve( numEntriesA+numEntriesB );
    const Int firstLocalRow = C.FirstLocalRow();
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueLocalUpdate( A.Row(e)-firstLocalRow, A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueLocalUpdate( B.Row(e)-firstLocalRow, B.Col(e)+nA, B.Value(e) );
    C.ProcessLocalQueues();
}

template<typename T>
void VCat
( const DistSparseMatrix<T>& A,
  const DistSparseMatrix<T>& B,
        DistSparseMatrix<T>& C )
{
    EL_DEBUG_CSE
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");
    /*
    if( A.Grid().Comm() != B.Grid().Comm() )
        LogicError("A and B had different communicators");
    */

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesB = B.NumLocalEntries();
    C.SetGrid( A.Grid() );
    C.Resize( mA+mB, n );
    Zero( C );
    C.Reserve( numEntriesA+numEntriesB, numEntriesA+numEntriesB );
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueUpdate( B.Row(e)+mA, B.Col(e), B.Value(e) );
    C.ProcessQueues();
}

template<typename T>
void HCat
( const DistMultiVec<T>& A,
  const DistMultiVec<T>& B,
        DistMultiVec<T>& C )
{
    EL_DEBUG_CSE
    if( A.Height() != B.Height() )
        LogicError("A and B must be the same height for HCat");

    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    C.SetGrid( A.Grid() );
    C.Resize( m, nA+nB );
    Zero( C );

    const Int localHeight = C.LocalHeight();
    const auto& ALoc = A.LockedMatrix();
    const auto& BLoc = B.LockedMatrix();
    auto& CLoc = C.Matrix();
    auto CLocL = CLoc( IR(0,localHeight), IR(0,nA)     );
    auto CLocR = CLoc( IR(0,localHeight), IR(nA,nA+nB) );
    CLocL = ALoc;
    CLocR = BLoc;
}

template<typename T>
void VCat
( const DistMultiVec<T>& A,
  const DistMultiVec<T>& B,
        DistMultiVec<T>& C )
{
    EL_DEBUG_CSE
    if( A.Width() != B.Width() )
        LogicError("A and B must be the same width for VCat");

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int mLocA = A.LocalHeight();
    const Int mLocB = B.LocalHeight();
    const Int n = A.Width();

    C.SetGrid( A.Grid() );
    C.Resize( mA+mB, n );
    Zero( C );
    C.Reserve( (mLocA+mLocB)*n );
    for( Int iLoc=0; iLoc<mLocA; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        for( Int j=0; j<n; ++j )
            C.QueueUpdate( i, j, A.GetLocal(iLoc,j) );
    }
    for( Int iLoc=0; iLoc<mLocB; ++iLoc )
    {
        const Int i = B.GlobalRow(iLoc) + mA;
        for( Int j=0; j<n; ++j )
            C.QueueUpdate( i, j, B.GetLocal(iLoc,j) );
    }
    C.ProcessQueues();
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void HCat \
  ( const Matrix<T>& A, \
    const Matrix<T>& B, \
          Matrix<T>& C ); \
  EL_EXTERN template void VCat \
  ( const Matrix<T>& A, \
    const Matrix<T>& B, \
          Matrix<T>& C ); \
  EL_EXTERN template void HCat \
  ( const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  EL_EXTERN template void VCat \
  ( const AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  EL_EXTERN template void HCat \
  ( const SparseMatrix<T>& A, \
    const SparseMatrix<T>& B, \
          SparseMatrix<T>& C ); \
  EL_EXTERN template void VCat \
  ( const SparseMatrix<T>& A, \
    const SparseMatrix<T>& B, \
          SparseMatrix<T>& C ); \
  EL_EXTERN template void HCat \
  ( const DistSparseMatrix<T>& A, \
    const DistSparseMatrix<T>& B, \
          DistSparseMatrix<T>& C ); \
  EL_EXTERN template void VCat \
  ( const DistSparseMatrix<T>& A, \
    const DistSparseMatrix<T>& B, \
          DistSparseMatrix<T>& C ); \
  EL_EXTERN template void HCat \
  ( const DistMultiVec<T>& A, \
    const DistMultiVec<T>& B, \
          DistMultiVec<T>& C ); \
  EL_EXTERN template void VCat \
  ( const DistMultiVec<T>& A, \
    const DistMultiVec<T>& B, \
          DistMultiVec<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_CONCATENATE_HPP
