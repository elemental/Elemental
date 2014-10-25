/*
   Copyright (c) 2009-2014, Jack Poulson
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
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A; 
}

template<typename S,typename T>
void Copy( const Matrix<S>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    auto convert = []( const S alpha ) { return T(alpha); };
    EntrywiseMap( A, B, std::function<T(S)>(convert) );
}

template<typename T,Dist U,Dist V>
inline void Copy( const AbstractDistMatrix<T>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
inline void Copy( const AbstractDistMatrix<S>& A, DistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
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
( const AbstractBlockDistMatrix<T>& A, BlockDistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
inline void Copy
( const AbstractBlockDistMatrix<S>& A, BlockDistMatrix<T,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
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
    BlockDistMatrix<S,U,V> BOrig(A.Grid());
    BOrig.AlignWith( B );
    BOrig = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BOrig.LockedMatrix(), B.Matrix() );
}

template<typename S,typename T>
void Copy( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

template<typename S,typename T>
void Copy( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define GUARD(CDIST,RDIST) B.ColDist() == CDIST && B.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& BCast = dynamic_cast<BlockDistMatrix<T,CDIST,RDIST>&>(B); \
        Copy( A, BCast );
    #include "El/macros/GuardAndPayload.h"
}

void Copy( const Graph& A, Graph& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy [Graph]"))
    B = A;
}

void Copy( const DistGraph& A, DistGraph& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy [DistGraph]"))
    B = A;
}

void CopyFromRoot( const DistGraph& distGraph, Graph& graph )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromRoot"))
    const mpi::Comm comm = distGraph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEdges = distGraph.NumLocalEdges();
    std::vector<int> edgeSizes(commSize), edgeOffsets(commSize);
    mpi::AllGather( &numLocalEdges, 1, &edgeSizes[0], 1, comm );
    int numEdges=0;
    for( int q=0; q<commSize; ++q )
    {
        edgeOffsets[q] = numEdges;
        numEdges += edgeSizes[q];
    }

    graph.Resize( distGraph.NumSources(), distGraph.NumTargets() );
    graph.Reserve( numEdges );
    graph.sources_.resize( numEdges );
    graph.targets_.resize( numEdges );
    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      graph.SourceBuffer(), &edgeSizes[0], &edgeOffsets[0], commRank, comm );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      graph.TargetBuffer(), &edgeSizes[0], &edgeOffsets[0], commRank, comm );
    graph.MakeConsistent();
}

void CopyFromNonRoot( const DistGraph& distGraph, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromRoot"))
    const mpi::Comm comm = distGraph.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEdges = distGraph.NumLocalEdges();
    std::vector<int> edgeSizes(commSize), edgeOffsets(commSize);
    mpi::AllGather( &numLocalEdges, 1, &edgeSizes[0], 1, comm );
    int numEdges=0;
    for( int q=0; q<commSize; ++q )
    {
        edgeOffsets[q] = numEdges;
        numEdges += edgeSizes[q];
    }

    mpi::Gather
    ( distGraph.LockedSourceBuffer(), numLocalEdges,
      (Int*)0, &edgeSizes[0], &edgeOffsets[0], root, comm );
    mpi::Gather
    ( distGraph.LockedTargetBuffer(), numLocalEdges,
      (Int*)0, &edgeSizes[0], &edgeOffsets[0], root, comm );
}

template<typename T>
void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy [SparseMatrix]"))
    B = A;
}

template<typename T>
void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy [DistSparseMatrix]"))
    B = A;
}

template<typename T>
void CopyFromRoot( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromRoot"))
    const mpi::Comm comm = ADist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEntries = ADist.NumLocalEntries();
    std::vector<int> entrySizes(commSize), entryOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &entrySizes[0], 1, comm );
    int numEntries=0;
    for( int q=0; q<commSize; ++q )
    {
        entryOffsets[q] = numEntries;
        numEntries += entrySizes[q];
    }

    A.Resize( ADist.Height(), ADist.Width() );
    A.Reserve( numEntries );
    A.graph_.sources_.resize( numEntries );
    A.graph_.targets_.resize( numEntries );
    A.vals_.resize( numEntries );
    mpi::Gather
    ( ADist.LockedSourceBuffer(), numLocalEntries,
      A.SourceBuffer(), &entrySizes[0], &entryOffsets[0], commRank, comm );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      A.TargetBuffer(), &entrySizes[0], &entryOffsets[0], commRank, comm );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      A.ValueBuffer(), &entrySizes[0], &entryOffsets[0], commRank, comm );
    A.MakeConsistent();
}

template<typename T>
void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromRoot"))
    const mpi::Comm comm = ADist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEntries = ADist.NumLocalEntries();
    std::vector<int> entrySizes(commSize), entryOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &entrySizes[0], 1, comm );
    int numEntries=0;
    for( int q=0; q<commSize; ++q )
    {
        entryOffsets[q] = numEntries;
        numEntries += entrySizes[q];
    }

    mpi::Gather
    ( ADist.LockedSourceBuffer(), numLocalEntries,
      (Int*)0, &entrySizes[0], &entryOffsets[0], root, comm );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      (Int*)0, &entrySizes[0], &entryOffsets[0], root, comm );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      (T*)0, &entrySizes[0], &entryOffsets[0], root, comm );
}

template<typename T>
void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy [DistMultiVec]"))
    B = A;
}

template<typename T>
void CopyFromRoot( const DistMultiVec<T>& XDist, Matrix<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromRoot"))
    const mpi::Comm comm = XDist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );

    const int numLocalEntries = XDist.LocalHeight()*XDist.Width();
    std::vector<int> entrySizes(commSize), entryOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &entrySizes[0], 1, comm );
    int numEntries=0;
    for( int q=0; q<commSize; ++q )
    {
        entryOffsets[q] = numEntries;
        numEntries += entrySizes[q];
    }

    X.Resize( XDist.Height(), XDist.Width(), XDist.Height() );
    const auto& XDistLoc = XDist.LockedMatrix();
    if( XDistLoc.Height() == XDistLoc.LDim() )
    {
        mpi::Gather
        ( XDistLoc.LockedBuffer(), numLocalEntries,
          X.Buffer(), &entrySizes[0], &entryOffsets[0], commRank, comm );
    }
    else
    {
        std::vector<T> sendBuf( numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc.Get(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          X.Buffer(), &entrySizes[0], &entryOffsets[0], commRank, comm );
    }
}

template<typename T>
void CopyFromNonRoot( const DistMultiVec<T>& XDist, Int root )
{
    DEBUG_ONLY(CallStackEntry cse("CopyFromNonRoot"))
    const mpi::Comm comm = XDist.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commRank == root )
        LogicError("Called CopyFromNonRoot from root");

    const int numLocalEntries = XDist.LocalHeight()*XDist.Width();
    std::vector<int> entrySizes(commSize), entryOffsets(commSize);
    mpi::AllGather( &numLocalEntries, 1, &entrySizes[0], 1, comm );
    int numEntries=0;
    for( int q=0; q<commSize; ++q )
    {
        entryOffsets[q] = numEntries;
        numEntries += entrySizes[q];
    }

    const auto& XDistLoc = XDist.LockedMatrix();
    if( XDistLoc.Height() == XDistLoc.LDim() )
    {
        mpi::Gather
        ( XDistLoc.LockedBuffer(), numLocalEntries,
          (T*)0, &entrySizes[0], &entryOffsets[0], root, comm );
    }
    else
    {
        std::vector<T> sendBuf( numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc.Get(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          (T*)0, &entrySizes[0], &entryOffsets[0], root, comm );
    }
}

// TODO: include guards so that certain datatypes can be properly disabled 
#define CONVERT(S,T) \
  template void Copy( const Matrix<S>& A, Matrix<T>& B ); \
  template void Copy \
  ( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B ); \
  template void Copy \
  ( const AbstractBlockDistMatrix<S>& A, AbstractBlockDistMatrix<T>& B );

#define SAME(T) \
  CONVERT(T,T) \
  template void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  template void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B ); \
  template void CopyFromRoot \
  ( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A ); \
  template void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, Int root ); \
  template void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B ); \
  template void CopyFromRoot( const DistMultiVec<T>& ADist, Matrix<T>& A ); \
  template void CopyFromNonRoot( const DistMultiVec<T>& ADist, Int root );

#define PROTO_INT(T) SAME(T) 

#define PROTO_REAL(Real) \
  SAME(Real) \
  /* Promotions up to Real */ \
  CONVERT(Int,Real) \
  /* Promotions up from Real */ \
  CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  SAME(C) \
  /* Promotions up to C */ \
  CONVERT(Int,C)

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  /* Promotions up from float */ \
  CONVERT(float,double) \
  CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  /* Promotions down to float */ \
  CONVERT(double,float) \
  /* Mixed conversion */ \
  CONVERT(double,Complex<float>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  /* Promotions up from Complex<float> */ \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  /* Promotions down from Complex<double> */ \
  CONVERT(Complex<double>,Complex<float>)

#include "El/macros/Instantiate.h"

} // namespace El
