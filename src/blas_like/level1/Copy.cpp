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
void Copy( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    const Int height = A.Height();
    const Int width = A.Width();
    B.Resize( height, width ); 

    lapack::Copy
    ( 'F', A.Height(), A.Width(),
      A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
}

template<typename S,typename T>
void Copy( const Matrix<S>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    EntrywiseMap( A, B, function<T(S)>(&Caster<S,T>::Cast) );
}

template<typename T,Dist U,Dist V>
inline void Copy( const ElementalMatrix<T>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
inline void Copy( const ElementalMatrix<S>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    if( A.Grid() == B.Grid() && A.ColDist() == U && A.RowDist() == V )
    {
        if( !B.RootConstrained() )
            B.SetRoot( A.Root() );
        if( !B.ColConstrained() )
            B.AlignCols( A.ColAlign() );
        if( !B.RowConstrained() )
            B.AlignRows( A.RowAlign() );
        if( A.Root() == B.Root() && 
            A.ColAlign() == B.ColAlign() && A.RowAlign() == B.RowAlign() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }
    DistMatrix<S,U,V> BOrig(A.Grid());
    BOrig.AlignWith( B );
    BOrig = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BOrig.LockedMatrix(), B.Matrix() );
}

template<typename T,Dist U,Dist V>
inline void Copy
( const BlockMatrix<T>& A, DistMatrix<T,U,V,BLOCK>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
inline void Copy
( const BlockMatrix<S>& A, DistMatrix<T,U,V,BLOCK>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    if( A.Grid() == B.Grid() && A.ColDist() == U && A.RowDist() == V )
    {
        if( !B.RootConstrained() )
            B.SetRoot( A.Root() );
        if( !B.ColConstrained() )
            B.AlignColsWith( A.DistData() );
        if( !B.RowConstrained() )
            B.AlignRowsWith( A.DistData() );
        if( A.Root() == B.Root() && 
            A.ColAlign() == B.ColAlign() && 
            A.RowAlign() == B.RowAlign() && 
            A.ColCut() == B.ColCut() &&
            A.RowCut() == B.RowCut() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }
    DistMatrix<S,U,V,BLOCK> BOrig(A.Grid());
    BOrig.AlignWith( B );
    BOrig = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BOrig.LockedMatrix(), B.Matrix() );
}

template<typename S,typename T>
void Copy( const ElementalMatrix<S>& A, ElementalMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename S,typename T>
void Copy( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    const DistWrap wrapA=A.Wrap(), wrapB=B.Wrap();
    if( wrapA == ELEMENT && wrapB == ELEMENT )
    {
        auto& ACast = dynamic_cast<const ElementalMatrix<T>&>(A);
        auto& BCast = dynamic_cast<ElementalMatrix<T>&>(B);
        BCast = ACast;
    }
    else if( wrapA == BLOCK && wrapB == BLOCK )
    {
        auto& ACast = dynamic_cast<const BlockMatrix<T>&>(A);
        auto& BCast = dynamic_cast<BlockMatrix<T>&>(B);
        BCast = ACast;
    }
    else 
    {
        Zeros( B, A.Height(), A.Width() );
        if( A.RedundantRank() == 0 )
        {
            const Int localHeight = A.LocalHeight();
            const Int localWidth = A.LocalWidth();
            B.Reserve( localHeight*localWidth );
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                {
                    const Int i = A.GlobalRow(iLoc);
                    B.QueueUpdate
                    ( i, j, Caster<S,T>::Cast(A.GetLocal(iLoc,jLoc)) );
                }
            }
        }
        const bool includeViewers = (A.Grid() != B.Grid());
        B.ProcessQueues( includeViewers );
    }
}

template<typename S,typename T>
void Copy( const BlockMatrix<S>& A, BlockMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
      auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST,BLOCK>&>(B); \
      Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

void Copy( const Graph& A, Graph& B )
{
    DEBUG_ONLY(CSE cse("Copy [Graph]"))
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();

    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.consistent_ = A.consistent_;
    B.sourceOffsets_ = A.sourceOffsets_;
    B.ProcessQueues();
}

void Copy( const Graph& A, DistGraph& B )
{
    DEBUG_ONLY(CSE cse("Copy [Graph/DistGraph]"))
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();

    B.SetComm( mpi::COMM_SELF );
    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.locallyConsistent_ = A.consistent_;
    B.localSourceOffsets_ = A.sourceOffsets_;
    B.ProcessLocalQueues();
}

void Copy( const DistGraph& A, Graph& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistGraph/Graph]"))
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();
    mpi::Comm comm = A.Comm();
    if( mpi::Size(comm) != 1 )
        LogicError("Cannot yet construct sequential graph from distributed");

    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.consistent_ = A.locallyConsistent_;
    B.sourceOffsets_ = A.localSourceOffsets_;
    B.ProcessQueues();
}

void Copy( const DistGraph& A, DistGraph& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistGraph]"))
    const Int numSources = A.NumSources();
    const Int numTargets = A.NumTargets();
    
    B.SetComm( A.Comm() );
    B.Resize( numSources, numTargets );
    // Directly assign instead of queueing up the individual edges
    B.sources_ = A.sources_;
    B.targets_ = A.targets_;
    B.locallyConsistent_ = A.locallyConsistent_;
    B.localSourceOffsets_ = A.localSourceOffsets_;
    B.ProcessLocalQueues();
}

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC>& B, bool includingViewers )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    if( B.CrossRank() != B.Root() )
        LogicError("Called CopyFromRoot from non-root");
    B.Resize( A.Height(), A.Width() );
    B.MakeSizeConsistent( includingViewers );
    B.Matrix() = A;
}

template<typename T>
void CopyFromNonRoot( DistMatrix<T,CIRC,CIRC>& B, bool includingViewers )
{
    DEBUG_ONLY(CSE cse("CopyFromNonRoot"))
    if( B.CrossRank() == B.Root() )
        LogicError("Called CopyFromNonRoot from root");
    B.MakeSizeConsistent( includingViewers );
}

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC,BLOCK>& B,
  bool includingViewers )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    if( B.CrossRank() != B.Root() )
        LogicError("Called CopyFromRoot from non-root");
    B.Resize( A.Height(), A.Width() );
    B.MakeSizeConsistent( includingViewers );
    B.Matrix() = A;
}

template<typename T>
void CopyFromNonRoot
( DistMatrix<T,CIRC,CIRC,BLOCK>& B, bool includingViewers )
{
    DEBUG_ONLY(CSE cse("CopyFromNonRoot"))
    if( B.CrossRank() == B.Root() )
        LogicError("Called CopyFromNonRoot from root");
    B.MakeSizeConsistent( includingViewers );
}

void CopyFromRoot( const DistGraph& distGraph, Graph& graph )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    const mpi::Comm comm = distGraph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEdges = distGraph.NumLocalEdges();
    vector<int> edgeSizes(commSize);
    mpi::AllGather( &numLocalEdges, 1, edgeSizes.data(), 1, comm );
    vector<int> edgeOffsets;
    const int numEdges = Scan( edgeSizes, edgeOffsets );

    graph.Resize( distGraph.NumSources(), distGraph.NumTargets() );
    graph.Reserve( numEdges );
    graph.sources_.resize( numEdges );
    graph.targets_.resize( numEdges );
    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      graph.SourceBuffer(), edgeSizes.data(), edgeOffsets.data(), 
      commRank, comm );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      graph.TargetBuffer(), edgeSizes.data(), edgeOffsets.data(), 
      commRank, comm );
    graph.ProcessQueues();
}

void CopyFromNonRoot( const DistGraph& distGraph, int root )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    const mpi::Comm comm = distGraph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEdges = distGraph.NumLocalEdges();
    vector<int> edgeSizes(commSize);
    mpi::AllGather( &numLocalEdges, 1, edgeSizes.data(), 1, comm );
    vector<int> edgeOffsets;
    Scan( edgeSizes, edgeOffsets );

    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      (Int*)0, edgeSizes.data(), edgeOffsets.data(), root, comm );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      (Int*)0, edgeSizes.data(), edgeOffsets.data(), root, comm );
}

template<typename T>
void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [SparseMatrix]"))
    B = A;
}

template<typename S,typename T>
void Copy( const SparseMatrix<S>& A, SparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    EntrywiseMap( A, B, function<T(S)>(&Caster<S,T>::Cast) );
}

template<typename S,typename T>
void Copy( const SparseMatrix<S>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    const S* AValBuf = A.LockedValueBuffer();
    const Int* ARowBuf = A.LockedSourceBuffer();
    const Int* AColBuf = A.LockedTargetBuffer();
    
    T* BBuf = B.Buffer();
    const Int BLDim = B.LDim();
    
    Zeros( B, m, n );
    for( Int e=0; e<numEntries; ++e )
        BBuf[ARowBuf[e]+AColBuf[e]*BLDim] = Caster<S,T>::Cast(AValBuf[e]);
}

template<typename T>
void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistSparseMatrix]"))
    B = A;
}

template<typename S,typename T>
void Copy( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    EntrywiseMap( A, B, function<T(S)>(&Caster<S,T>::Cast) );
}

template<typename S,typename T>
void Copy( const DistSparseMatrix<S>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumLocalEntries();
    Zeros( B, m, n );
    B.Reserve( numEntries );
    for( Int e=0; e<numEntries; ++e )
        B.QueueUpdate( A.Row(e), A.Col(e), Caster<S,T>::Cast(A.Value(e)) );
    B.ProcessQueues();
}

template<typename T>
void CopyFromRoot( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    const mpi::Comm comm = ADist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEntries = ADist.NumLocalEntries();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, comm );
    vector<int> entryOffs;
    const int numEntries = Scan( entrySizes, entryOffs );

    A.Resize( ADist.Height(), ADist.Width() );
    A.Reserve( numEntries );
    A.graph_.sources_.resize( numEntries );
    A.graph_.targets_.resize( numEntries );
    A.vals_.resize( numEntries );
    mpi::Gather
    ( ADist.LockedSourceBuffer(), numLocalEntries,
      A.SourceBuffer(), entrySizes.data(), entryOffs.data(), 
      commRank, comm );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      A.TargetBuffer(), entrySizes.data(), entryOffs.data(), 
      commRank, comm );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      A.ValueBuffer(), entrySizes.data(), entryOffs.data(), 
      commRank, comm );
    A.ProcessQueues();
}

template<typename T>
void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, int root )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    const mpi::Comm comm = ADist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEntries = ADist.NumLocalEntries();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, comm );
    vector<int> entryOffs;
    Scan( entrySizes, entryOffs );

    mpi::Gather
    ( ADist.LockedSourceBuffer(), numLocalEntries,
      (Int*)0, entrySizes.data(), entryOffs.data(), root, comm );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      (Int*)0, entrySizes.data(), entryOffs.data(), root, comm );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      (T*)0, entrySizes.data(), entryOffs.data(), root, comm );
}

template<typename T>
void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistMultiVec]"))
    B.SetComm( A.Comm() );
    B.Resize( A.Height(), A.Width() );
    B.Matrix() = A.LockedMatrix();
}

template<typename S,typename T>
void Copy( const DistMultiVec<S>& A, DistMultiVec<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistMultiVec]"))
    EntrywiseMap( A, B, function<T(S)>(&Caster<S,T>::Cast) );
}

template<typename T>
void Copy( const DistMultiVec<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [DistMultiVec -> ADM]"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    Zeros( B, m, n );
    B.Reserve( mLoc*n );
    for( Int iLoc=0; iLoc<mLoc; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        for( Int j=0; j<n; ++j )
            B.QueueUpdate( i, j, A.GetLocal(iLoc,j) );
    }
    B.ProcessQueues();
}

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, DistMultiVec<T>& B )
{
    DEBUG_ONLY(CSE cse("Copy [ADM -> DistMultiVec]"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    mpi::Comm comm = A.Grid().Comm();
    B.SetComm( comm );
    Zeros( B, m, n );
    B.Reserve( mLoc*nLoc );
    for( Int iLoc=0; iLoc<mLoc; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            B.QueueUpdate( i, j, A.GetLocal(iLoc,jLoc) );
        }
    }
    B.ProcessQueues();
}

template<typename T>
void CopyFromRoot( const DistMultiVec<T>& XDist, Matrix<T>& X )
{
    DEBUG_ONLY(CSE cse("CopyFromRoot"))
    const Int m = XDist.Height();
    const Int n = XDist.Width();
    X.Resize( m, n, Max(m,1) );
    if( Min(m,n) == 0 )
        return;

    const mpi::Comm comm = XDist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEntries = XDist.LocalHeight()*n;
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, comm );
    vector<int> entryOffs;
    const int numEntries = Scan( entrySizes, entryOffs );

    vector<T> recvBuf;
    FastResize( recvBuf, numEntries );

    const auto& XDistLoc = XDist.LockedMatrix();
    if( XDistLoc.Height() == XDistLoc.LDim() )
    {
        mpi::Gather
        ( XDistLoc.LockedBuffer(), numLocalEntries,
          recvBuf.data(), entrySizes.data(), entryOffs.data(), 
          commRank, comm );
    }
    else
    {
        vector<T> sendBuf;
        FastResize( sendBuf, numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc.Get(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          recvBuf.data(), entrySizes.data(), entryOffs.data(), 
          commRank, comm );
    }
    for( Int q=0; q<commSize; ++q )
    {
        const Int iOff = entryOffs[q]/n;
        const Int iSize = entrySizes[q]/n;
        for( Int t=0; t<entrySizes[q]; ++t )
            X.Set( iOff+(t%iSize), t/iSize, recvBuf[entryOffs[q]+t] );
    }
}

template<typename T>
void CopyFromNonRoot( const DistMultiVec<T>& XDist, int root )
{
    DEBUG_ONLY(CSE cse("CopyFromNonRoot"))
    const Int m = XDist.Height();
    const Int n = XDist.Width();
    if( Min(m,n) == 0 )
        return;

    const mpi::Comm comm = XDist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Called CopyFromNonRoot from root");

    const int numLocalEntries = XDist.LocalHeight()*XDist.Width();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, comm );
    vector<int> entryOffs;
    Scan( entrySizes, entryOffs );

    const auto& XDistLoc = XDist.LockedMatrix();
    if( XDistLoc.Height() == XDistLoc.LDim() )
    {
        mpi::Gather
        ( XDistLoc.LockedBuffer(), numLocalEntries,
          (T*)0, entrySizes.data(), entryOffs.data(), root, comm );
    }
    else
    {
        vector<T> sendBuf;
        FastResize( sendBuf, numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc.Get(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          (T*)0, entrySizes.data(), entryOffs.data(), root, comm );
    }
}

#define CONVERT(S,T) \
  template void Copy( const Matrix<S>& A, Matrix<T>& B ); \
  template void Copy \
  ( const DistSparseMatrix<S>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const ElementalMatrix<S>& A, ElementalMatrix<T>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const BlockMatrix<S>& A, BlockMatrix<T>& B ); \
  template void Copy( const SparseMatrix<S>& A, SparseMatrix<T>& B ); \
  template void Copy( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B ); \
  template void Copy( const SparseMatrix<S>& A, Matrix<T>& B ); \
  template void Copy( const DistMultiVec<S>& A, DistMultiVec<T>& B );

#define SAME(T) \
  CONVERT(T,T) \
  template void CopyFromRoot \
  ( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC>& B, bool includingViewers ); \
  template void CopyFromNonRoot \
  ( DistMatrix<T,CIRC,CIRC>& B, bool includingViewers ); \
  template void CopyFromRoot \
  ( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A ); \
  template void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, int root ); \
  template void Copy( const DistMultiVec<T>& A, AbstractDistMatrix<T>& B ); \
  template void Copy( const AbstractDistMatrix<T>& A, DistMultiVec<T>& B ); \
  template void CopyFromRoot( const DistMultiVec<T>& ADist, Matrix<T>& A ); \
  template void CopyFromNonRoot( const DistMultiVec<T>& ADist, int root );

#define PROTO_INT(T) SAME(T) 

#define PROTO_REAL(Real) \
  SAME(Real) \
  CONVERT(Int,Real) \
  CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  SAME(C) \
  CONVERT(Int,C)

#ifdef EL_HAVE_QUAD

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Quad) \
  CONVERT(float,Complex<double>) \
  CONVERT(float,Complex<Quad>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Quad) \
  CONVERT(double,Complex<float>) \
  CONVERT(double,Complex<Quad>)

#define PROTO_QUAD \
  PROTO_REAL(Quad) \
  CONVERT(Quad,float) \
  CONVERT(Quad,double) \
  CONVERT(Quad,Complex<float>) \
  CONVERT(Quad,Complex<double>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>) \
  CONVERT(Complex<float>,Complex<Quad>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>) \
  CONVERT(Complex<double>,Complex<Quad>)

#define PROTO_COMPLEX_QUAD \
  PROTO_COMPLEX(Complex<Quad>) \
  CONVERT(Complex<Quad>,Complex<float>) \
  CONVERT(Complex<Quad>,Complex<double>)

#else

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Complex<float>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>)

#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
