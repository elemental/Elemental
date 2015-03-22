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
void HCat( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("HCat"))
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");
    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    Zeros( C, m, nA+nB );
    auto CL = C( IR(0,m), IR(0,nA)     );
    auto CR = C( IR(0,m), IR(nA,nA+nB) );
    CL = A;
    CR = B;
}

template<typename T>
void VCat( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("VCat"))
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    Zeros( C, mA+mB, n );
    auto CT = C( IR(0,mA),     IR(0,n) );
    auto CB = C( IR(mA,mA+mB), IR(0,n) );
    CT = A;
    CB = B;
}

template<typename T>
inline void HCat
( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, 
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat");
    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();

    auto CPtr = WriteProxy<T,MC,MR>(&CPre);
    auto& C = *CPtr;

    Zeros( C, m, nA+nB );
    auto CL = C( IR(0,m), IR(0,nA)     );
    auto CR = C( IR(0,m), IR(nA,nA+nB) );
    CL = A;
    CR = B;
}

template<typename T>
void VCat
( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, 
        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(CallStackEntry cse("VCat"))
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat");
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();

    auto CPtr = WriteProxy<T,MC,MR>(&CPre);
    auto& C = *CPtr;

    Zeros( C, mA+mB, n );
    auto CT = C( IR(0,mA),     IR(0,n) );
    auto CB = C( IR(mA,mA+mB), IR(0,n) );
    CT = A;
    CB = B;
}

template<typename T>
void HCat
( const SparseMatrix<T>& A, const SparseMatrix<T>& B, 
        SparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("HCat"))
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat"); 

    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();
    
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    Zeros( C, m, nA+nB );
    C.Reserve( numEntriesA+numEntriesB ); 
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueUpdate( B.Row(e), B.Col(e)+nA, B.Value(e) );
    C.MakeConsistent();
}

template<typename T>
void VCat
( const SparseMatrix<T>& A, const SparseMatrix<T>& B, 
        SparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("VCat"))
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat"); 

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesB = B.NumEntries();
    Zeros( C, mA+mB, n );
    C.Reserve( numEntriesA+numEntriesB ); 
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueUpdate( A.Row(e), A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueUpdate( B.Row(e)+mA, B.Col(e), B.Value(e) );
    C.MakeConsistent();
}

template<typename T>
void HCat
( const DistSparseMatrix<T>& A, const DistSparseMatrix<T>& B, 
        DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("HCat"))
    if( A.Height() != B.Height() )
        LogicError("Incompatible heights for HCat"); 
    /*
    if( A.Comm() != B.Comm() )
        LogicError("A and B had different communicators");
    */

    const Int m = A.Height();
    const Int nA = A.Width();
    const Int nB = B.Width();
    
    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesB = B.NumLocalEntries();
    C.SetComm( A.Comm() );
    Zeros( C, m, nA+nB );
    C.Reserve( numEntriesA+numEntriesB ); 
    const Int firstLocalRow = C.FirstLocalRow();
    for( Int e=0; e<numEntriesA; ++e )
        C.QueueLocalUpdate( A.Row(e)-firstLocalRow, A.Col(e), A.Value(e) );
    for( Int e=0; e<numEntriesB; ++e )
        C.QueueLocalUpdate( B.Row(e)-firstLocalRow, B.Col(e)+nA, B.Value(e) );
    C.MakeConsistent();
}

template<typename T>
void VCat
( const DistSparseMatrix<T>& A, const DistSparseMatrix<T>& B, 
        DistSparseMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("VCat"))
    if( A.Width() != B.Width() )
        LogicError("Incompatible widths for VCat"); 
    /*
    if( A.Comm() != B.Comm() )
        LogicError("A and B had different communicators");
    */

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    
    const Int numEntriesA = A.NumLocalEntries();
    const Int numEntriesB = B.NumLocalEntries();
    C.SetComm( A.Comm() );
    Zeros( C, mA+mB, n );

    // Compute the metadata
    // --------------------
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEntriesA; ++e )
        ++sendCounts[ C.RowOwner(A.Row(e)) ];
    for( Int e=0; e<numEntriesB; ++e )
        ++sendCounts[ C.RowOwner(B.Row(e)+mA) ];
    vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );
    vector<int> sendOffs, recvOffs; 
    const int totalSend = Scan( sendCounts, sendOffs );
    const int totalRecv = Scan( recvCounts, recvOffs );
    // Pack
    // ----
    vector<ValueIntPair<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int e=0; e<numEntriesA; ++e )
    {
        const Int i = A.Row(e);    
        const Int j = A.Col(e);
        const T value = A.Value(e);
        const int owner = C.RowOwner(i);
        sendBuf[offs[owner]].value = value;
        sendBuf[offs[owner]].indices[0] = i;
        sendBuf[offs[owner]].indices[1] = j;
        ++offs[owner];
    }
    for( Int e=0; e<numEntriesB; ++e )
    {
        const Int i = B.Row(e)+mA;
        const Int j = B.Col(e);
        const T value = B.Value(e);
        const int owner = C.RowOwner(i);
        sendBuf[offs[owner]].value = value;
        sendBuf[offs[owner]].indices[0] = i;
        sendBuf[offs[owner]].indices[1] = j;
        ++offs[owner];
    }
    // Exchange
    // --------
    vector<ValueIntPair<T>> recvBuf(totalRecv);
    mpi::AllToAll
    ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
      recvBuf.data(), recvCounts.data(), recvOffs.data(), comm );
    // Unpack
    // ------
    C.Reserve( totalRecv );
    for( Int e=0; e<totalRecv; ++e )
        C.QueueLocalUpdate
        ( recvBuf[e].indices[0]-C.FirstLocalRow(), recvBuf[e].indices[1], 
          recvBuf[e].value );
    C.MakeConsistent();
}

#define PROTO(T) \
  template void HCat \
  ( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ); \
  template void VCat \
  ( const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C ); \
  template void HCat \
  ( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  template void VCat \
  ( const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
          AbstractDistMatrix<T>& C ); \
  template void HCat \
  ( const SparseMatrix<T>& A, const SparseMatrix<T>& B, SparseMatrix<T>& C ); \
  template void VCat \
  ( const SparseMatrix<T>& A, const SparseMatrix<T>& B, SparseMatrix<T>& C ); \
  template void HCat \
  ( const DistSparseMatrix<T>& A, const DistSparseMatrix<T>& B, \
          DistSparseMatrix<T>& C ); \
  template void VCat \
  ( const DistSparseMatrix<T>& A, const DistSparseMatrix<T>& B, \
          DistSparseMatrix<T>& C );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
