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
  const vector<Int>& I, const vector<Int>& J, 
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
  const vector<Int>& I, const vector<Int>& J, 
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
  const vector<Int>& I, const vector<Int>& J, 
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
( const Matrix<T>& A, const vector<Int>& I, const vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    Matrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>> GetRealPartOfSubmatrix
( const Matrix<T>& A, const vector<Int>& I, const vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    Matrix<Base<T>> ASub;
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
Matrix<Base<T>> GetImagPartOfSubmatrix
( const Matrix<T>& A, const vector<Int>& I, const vector<Int>& J )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    Matrix<Base<T>> ASub;
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
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
  const vector<Int>& I, const vector<Int>& J, 
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
  const vector<Int>& I, const vector<Int>& J, 
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
  const vector<Int>& I, const vector<Int>& J )
{
    DistMatrix<T,STAR,STAR> ASub( A.Grid() );
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J )
{
    DistMatrix<Base<T>,STAR,STAR> ASub( A.Grid() );
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J )
{
    DistMatrix<Base<T>,STAR,STAR> ASub( A.Grid() );
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        SparseMatrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    Zeros( ASub, mSub, nSub );

    // Reserve the number of nonzeros that live within the submatrix
    Int numNonzerosSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ++numNonzerosSub;
        }
    }
    ASub.Reserve( numNonzerosSub );

    // Insert the nonzeros
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ASub.QueueUpdate( i-I.beg, j-J.beg, A.Value(e) );
        }
    }
    ASub.MakeConsistent();
}

template<typename T>
void GetRealPartOfSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        SparseMatrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    Zeros( ASub, mSub, nSub );

    // Reserve the number of nonzeros that live within the submatrix
    Int numNonzerosSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ++numNonzerosSub;
        }
    }
    ASub.Reserve( numNonzerosSub );

    // Insert the nonzeros
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ASub.QueueUpdate( i-I.beg, j-J.beg, RealPart(A.Value(e)) );
        }
    }
    ASub.MakeConsistent();
}

template<typename T>
void GetImagPartOfSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        SparseMatrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    Zeros( ASub, mSub, nSub );

    // Reserve the number of nonzeros that live within the submatrix
    Int numNonzerosSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ++numNonzerosSub;
        }
    }
    ASub.Reserve( numNonzerosSub );

    // Insert the nonzeros
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int rowOff = A.EntryOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ASub.QueueUpdate( i-I.beg, j-J.beg, ImagPart(A.Value(e)) );
        }
    }
    ASub.MakeConsistent();
}

template<typename T>
SparseMatrix<T> GetSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    SparseMatrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
SparseMatrix<Base<T>> GetRealPartOfSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    SparseMatrix<Base<T>> ASub;
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
SparseMatrix<Base<T>> GetImagPartOfSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    SparseMatrix<Base<T>> ASub;
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        DistSparseMatrix<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int numEntries = A.NumLocalEntries();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++sendCounts[ ASub.RowOwner(i-I.beg) ];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<T>> sendBuf(totalSend);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const T value = A.Value(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
        {
            const Int iSub = i - I.beg;
            const Int jSub = j - J.beg;
            const int owner = ASub.RowOwner( iSub );
            sendBuf[offs[owner]++] = Entry<T>{ iSub, jSub, value };
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    ASub.Reserve( recvBuf.size() );
    for( auto& entry : recvBuf )
        ASub.QueueUpdate( entry );
    ASub.MakeConsistent();
}

template<typename T>
void GetRealPartOfSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        DistSparseMatrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int numEntries = A.NumLocalEntries();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++sendCounts[ ASub.RowOwner(i-I.beg) ];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<Base<T>>> sendBuf(totalSend);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Base<T> value = RealPart(A.Value(e));
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
        {
            const Int iSub = i - I.beg;
            const Int jSub = j - J.beg;
            const int owner = ASub.RowOwner( iSub );
            sendBuf[offs[owner]++] = Entry<Base<T>>{ iSub, jSub, value };
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    ASub.Reserve( recvBuf.size() );
    for( auto& entry : recvBuf )
        ASub.QueueUpdate( entry );
    ASub.MakeConsistent();
}

template<typename T>
void GetImagPartOfSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        DistSparseMatrix<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int numEntries = A.NumLocalEntries();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++sendCounts[ ASub.RowOwner(i-I.beg) ];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<Base<T>>> sendBuf(totalSend);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Base<T> value = ImagPart(A.Value(e));
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
        {
            const Int iSub = i - I.beg;
            const Int jSub = j - J.beg;
            const int owner = ASub.RowOwner( iSub );
            sendBuf[offs[owner]++] = Entry<Base<T>>{ iSub, jSub, value };
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    ASub.Reserve( recvBuf.size() );
    for( auto& entry : recvBuf )
        ASub.QueueUpdate( entry );
    ASub.MakeConsistent();
}

template<typename T>
DistSparseMatrix<T> GetSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    DistSparseMatrix<T> ASub(A.Comm());
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistSparseMatrix<Base<T>> GetRealPartOfSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    DistSparseMatrix<Base<T>> ASub(A.Comm());
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistSparseMatrix<Base<T>> GetImagPartOfSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    DistSparseMatrix<Base<T>> ASub(A.Comm());
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J,
        DistMultiVec<T>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int localHeight = A.LocalHeight();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );
    
    // If no communication is necessary, take the easy and fast approach
    if( mSub == A.Height() )
    {
        for( Int j=J.beg; j<J.end; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ASub.SetLocal( iLoc, j-J.beg, A.GetLocal(iLoc,j) );
        return;
    }

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
            sendCounts[ ASub.RowOwner(i-I.beg) ] += J.end-J.beg;
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<T>> sendBuf(totalSend);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
        {
            for( Int j=J.beg; j<J.end; ++j )
            {
                const T value = A.GetLocal(iLoc,j);
                const Int iSub = i - I.beg;
                const Int jSub = j - J.beg;
                const int owner = ASub.RowOwner( iSub );
                sendBuf[offs[owner]++] = Entry<T>{ iSub, jSub, value };
            }
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    for( auto& entry : recvBuf )
        ASub.Set( entry );
}

template<typename T>
void GetRealPartOfSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J,
        DistMultiVec<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetRealPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int localHeight = A.LocalHeight();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
            sendCounts[ ASub.RowOwner(i-I.beg) ] += J.end-J.beg;
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<Base<T>>> sendBuf(totalSend);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
        {
            for( Int j=J.beg; j<J.end; ++j )
            {
                const Base<T> value = RealPart(A.GetLocal(iLoc,j));
                const Int iSub = i - I.beg;
                const Int jSub = j - J.beg;
                const int owner = ASub.RowOwner( iSub );
                sendBuf[offs[owner]++] = Entry<Base<T>>{ iSub, jSub, value };
            }
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    for( auto& entry : recvBuf )
        ASub.Set( entry );
}

template<typename T>
void GetImagPartOfSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J,
        DistMultiVec<Base<T>>& ASub )
{
    DEBUG_ONLY(CallStackEntry cse("GetImagPartOfSubmatrix"))
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int localHeight = A.LocalHeight();

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    ASub.SetComm( comm );
    Zeros( ASub, mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
            sendCounts[ ASub.RowOwner(i-I.beg) ] += J.end-J.beg;
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<Base<T>>> sendBuf(totalSend);
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
        {
            for( Int j=J.beg; j<J.end; ++j )
            {
                const Base<T> value = ImagPart(A.GetLocal(iLoc,j));
                const Int iSub = i - I.beg;
                const Int jSub = j - J.beg;
                const int owner = ASub.RowOwner( iSub );
                sendBuf[offs[owner]++] = Entry<Base<T>>{ iSub, jSub, value };
            }
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    for( auto& entry : recvBuf )
        ASub.Set( entry );
}

template<typename T>
DistMultiVec<T> GetSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J )
{
    DistMultiVec<T> ASub(A.Comm());
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMultiVec<Base<T>> GetRealPartOfSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J )
{
    DistMultiVec<Base<T>> ASub(A.Comm());
    GetRealPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
DistMultiVec<Base<T>> GetImagPartOfSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J )
{
    DistMultiVec<Base<T>> ASub(A.Comm());
    GetImagPartOfSubmatrix( A, I, J, ASub );
    return ASub;
}

#define PROTO(T) \
  template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          Matrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          Matrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          Matrix<Base<T>>& ASub ); \
  template Matrix<T> GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template Matrix<Base<T>> GetRealPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template Matrix<Base<T>> GetImagPartOfSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          AbstractDistMatrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          AbstractDistMatrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          AbstractDistMatrix<Base<T>>& ASub ); \
  template DistMatrix<T,STAR,STAR> GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template void GetSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          SparseMatrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          SparseMatrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          SparseMatrix<Base<T>>& ASub ); \
  template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          DistSparseMatrix<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          DistSparseMatrix<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          DistSparseMatrix<Base<T>>& ASub ); \
  template SparseMatrix<T> GetSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template SparseMatrix<Base<T>> GetRealPartOfSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template SparseMatrix<Base<T>> GetImagPartOfSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template DistSparseMatrix<T> GetSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template DistSparseMatrix<Base<T>> GetRealPartOfSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template DistSparseMatrix<Base<T>> GetImagPartOfSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template void GetSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J, \
          DistMultiVec<T>& ASub ); \
  template void GetRealPartOfSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J, \
          DistMultiVec<Base<T>>& ASub ); \
  template void GetImagPartOfSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J, \
          DistMultiVec<Base<T>>& ASub ); \
  template DistMultiVec<T> GetSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J ); \
  template DistMultiVec<Base<T>> GetRealPartOfSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J ); \
  template DistMultiVec<Base<T>> GetImagPartOfSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
