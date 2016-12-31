/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_HPP
#define EL_BLAS_COPY_HPP

#include <El/blas_like/level1/Copy/internal_decl.hpp>
#include <El/blas_like/level1/Copy/GeneralPurpose.hpp>
#include <El/blas_like/level1/Copy/util.hpp>

namespace El {

template<typename T>
void Copy( const Matrix<T>& A, Matrix<T>& B )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    B.Resize( height, width );

    lapack::Copy
    ( 'F', A.Height(), A.Width(),
      A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const Matrix<S>& A, Matrix<T>& B )
{
    EL_DEBUG_CSE
    EntrywiseMap( A, B, MakeFunction(Caster<S,T>::Cast) );
}

template<typename T,Dist U,Dist V>
void Copy( const ElementalMatrix<T>& A, DistMatrix<T,U,V>& B )
{
    EL_DEBUG_CSE
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
void Copy( const ElementalMatrix<S>& A, DistMatrix<T,U,V>& B )
{
    EL_DEBUG_CSE
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
void Copy( const BlockMatrix<T>& A, DistMatrix<T,U,V,BLOCK>& B )
{
    EL_DEBUG_CSE
    B = A;
}

// Datatype conversions should not be very common, and so it is likely best to
// avoid explicitly instantiating every combination
template<typename S,typename T,Dist U,Dist V>
void Copy( const BlockMatrix<S>& A, DistMatrix<T,U,V,BLOCK>& B )
{
    EL_DEBUG_CSE
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

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const ElementalMatrix<S>& A, ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    #define GUARD(CDIST,RDIST,WRAP) \
      B.ColDist() == CDIST && B.RowDist() == RDIST && ELEMENT == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
        auto& BCast = static_cast<DistMatrix<T,CDIST,RDIST,ELEMENT>&>(B); \
        Copy( A, BCast );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    const DistWrap wrapA=A.Wrap(), wrapB=B.Wrap();
    if( wrapA == ELEMENT && wrapB == ELEMENT )
    {
        auto& ACast = static_cast<const ElementalMatrix<T>&>(A);
        auto& BCast = static_cast<ElementalMatrix<T>&>(B);
        Copy( ACast, BCast );
    }
    else if( wrapA == BLOCK && wrapB == BLOCK )
    {
        auto& ACast = static_cast<const BlockMatrix<T>&>(A);
        auto& BCast = static_cast<BlockMatrix<T>&>(B);
        Copy( ACast, BCast );
    }
    else
    {
        copy::GeneralPurpose( A, B );
    }
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const AbstractDistMatrix<S>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    const DistWrap wrapA=A.Wrap(), wrapB=B.Wrap();
    if( wrapA == ELEMENT && wrapB == ELEMENT )
    {
        auto& ACast = static_cast<const ElementalMatrix<S>&>(A);
        auto& BCast = static_cast<ElementalMatrix<T>&>(B);
        Copy( ACast, BCast );
    }
    else if( wrapA == BLOCK && wrapB == BLOCK )
    {
        auto& ACast = static_cast<const BlockMatrix<S>&>(A);
        auto& BCast = static_cast<BlockMatrix<T>&>(B);
        Copy( ACast, BCast );
    }
    else
    {
        copy::GeneralPurpose( A, B );
    }
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const BlockMatrix<S>& A, BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    #define GUARD(CDIST,RDIST,WRAP) \
      B.ColDist() == CDIST && B.RowDist() == RDIST && BLOCK == WRAP
    #define PAYLOAD(CDIST,RDIST,WRAP) \
      auto& BCast = static_cast<DistMatrix<T,CDIST,RDIST,BLOCK>&>(B); \
      Copy( A, BCast );
    #include <El/macros/GuardAndPayload.h>
}

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC>& B, bool includingViewers )
{
    EL_DEBUG_CSE
    if( B.CrossRank() != B.Root() )
        LogicError("Called CopyFromRoot from non-root");
    B.Resize( A.Height(), A.Width() );
    B.MakeSizeConsistent( includingViewers );
    B.Matrix() = A;
}

template<typename T>
void CopyFromNonRoot( DistMatrix<T,CIRC,CIRC>& B, bool includingViewers )
{
    EL_DEBUG_CSE
    if( B.CrossRank() == B.Root() )
        LogicError("Called CopyFromNonRoot from root");
    B.MakeSizeConsistent( includingViewers );
}

template<typename T>
void CopyFromRoot
( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC,BLOCK>& B,
  bool includingViewers )
{
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    if( B.CrossRank() == B.Root() )
        LogicError("Called CopyFromNonRoot from root");
    B.MakeSizeConsistent( includingViewers );
}

template<typename T>
void Copy( const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    B = A;
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const SparseMatrix<S>& A, SparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    EntrywiseMap( A, B, MakeFunction(Caster<S,T>::Cast) );
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const SparseMatrix<S>& A, Matrix<T>& B )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    const S* AValBuf = A.LockedValueBuffer();
    const Int* ARowBuf = A.LockedSourceBuffer();
    const Int* AColBuf = A.LockedTargetBuffer();

    T* BBuf = B.Buffer();
    const Int BLDim = B.LDim();

    B.Resize( m, n );
    Zero( B );
    for( Int e=0; e<numEntries; ++e )
        BBuf[ARowBuf[e]+AColBuf[e]*BLDim] = Caster<S,T>::Cast(AValBuf[e]);
}

template<typename T>
void Copy( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    B = A;
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const DistSparseMatrix<S>& A, DistSparseMatrix<T>& B )
{
    EL_DEBUG_CSE
    EntrywiseMap( A, B, MakeFunction(Caster<S,T>::Cast) );
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const DistSparseMatrix<S>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumLocalEntries();
    B.Resize( m, n );
    Zero( B );
    B.Reserve( numEntries );
    for( Int e=0; e<numEntries; ++e )
        B.QueueUpdate( A.Row(e), A.Col(e), Caster<S,T>::Cast(A.Value(e)) );
    B.ProcessQueues();
}

template<typename T>
void CopyFromRoot( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A )
{
    EL_DEBUG_CSE
    const Grid& grid = ADist.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();

    const int numLocalEntries = ADist.NumLocalEntries();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, grid.Comm() );
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
      commRank, grid.Comm() );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      A.TargetBuffer(), entrySizes.data(), entryOffs.data(),
      commRank, grid.Comm() );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      A.ValueBuffer(), entrySizes.data(), entryOffs.data(),
      commRank, grid.Comm() );
    A.ProcessQueues();
}

template<typename T>
void CopyFromNonRoot( const DistSparseMatrix<T>& ADist, int root )
{
    EL_DEBUG_CSE
    const Grid& grid = ADist.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();
    if( commRank == root )
        LogicError("Root called CopyFromNonRoot");

    const int numLocalEntries = ADist.NumLocalEntries();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, grid.Comm() );
    vector<int> entryOffs;
    Scan( entrySizes, entryOffs );

    mpi::Gather
    ( ADist.LockedSourceBuffer(), numLocalEntries,
      (Int*)0, entrySizes.data(), entryOffs.data(), root, grid.Comm() );
    mpi::Gather
    ( ADist.LockedTargetBuffer(), numLocalEntries,
      (Int*)0, entrySizes.data(), entryOffs.data(), root, grid.Comm() );
    mpi::Gather
    ( ADist.LockedValueBuffer(), numLocalEntries,
      (T*)0, entrySizes.data(), entryOffs.data(), root, grid.Comm() );
}

template<typename T>
void Copy( const DistMultiVec<T>& A, DistMultiVec<T>& B )
{
    EL_DEBUG_CSE
    B.SetGrid( A.Grid() );
    B.Resize( A.Height(), A.Width() );
    B.Matrix() = A.LockedMatrix();
}

template<typename S,typename T,
         typename/*=EnableIf<CanCast<S,T>>*/>
void Copy( const DistMultiVec<S>& A, DistMultiVec<T>& B )
{
    EL_DEBUG_CSE
    EntrywiseMap( A, B, MakeFunction(Caster<S,T>::Cast) );
}

template<typename T>
void Copy( const DistMultiVec<T>& A, AbstractDistMatrix<T>& B )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    B.Resize( m, n );
    Zero( B );
    B.Reserve( mLoc*n );
    auto& ALoc = A.LockedMatrix();
    for( Int iLoc=0; iLoc<mLoc; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        for( Int j=0; j<n; ++j )
            B.QueueUpdate( i, j, ALoc(iLoc,j) );
    }
    B.ProcessQueues();
}

template<typename T>
void Copy( const AbstractDistMatrix<T>& A, DistMultiVec<T>& B )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    B.SetGrid( A.Grid() );
    B.Resize( m, n );
    Zero( B );
    B.Reserve( mLoc*nLoc );
    auto& ALoc = A.LockedMatrix();
    for( Int iLoc=0; iLoc<mLoc; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        for( Int jLoc=0; jLoc<nLoc; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            B.QueueUpdate( i, j, ALoc(iLoc,jLoc) );
        }
    }
    B.ProcessQueues();
}

template<typename T>
void CopyFromRoot( const DistMultiVec<T>& XDist, Matrix<T>& X )
{
    EL_DEBUG_CSE
    const Int m = XDist.Height();
    const Int n = XDist.Width();
    X.Resize( m, n, Max(m,1) );
    if( Min(m,n) == 0 )
        return;

    const Grid& grid = XDist.Grid();
    Output("grid.Size()=",grid.Size());
    const int commSize = grid.Size();
    const int commRank = grid.Rank();

    const int numLocalEntries = XDist.LocalHeight()*n;
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, grid.Comm() );
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
          commRank, grid.Comm() );
    }
    else
    {
        vector<T> sendBuf;
        FastResize( sendBuf, numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          recvBuf.data(), entrySizes.data(), entryOffs.data(),
          commRank, grid.Comm() );
    }
    for( Int q=0; q<commSize; ++q )
    {
        const Int iOff = entryOffs[q]/n;
        const Int iSize = entrySizes[q]/n;
        for( Int t=0; t<entrySizes[q]; ++t )
            X( iOff+(t%iSize), t/iSize ) = recvBuf[entryOffs[q]+t];
    }
}

template<typename T>
void CopyFromNonRoot( const DistMultiVec<T>& XDist, int root )
{
    EL_DEBUG_CSE
    const Int m = XDist.Height();
    const Int n = XDist.Width();
    if( Min(m,n) == 0 )
        return;

    const Grid& grid = XDist.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();
    if( commRank == root )
        LogicError("Called CopyFromNonRoot from root");

    const int numLocalEntries = XDist.LocalHeight()*XDist.Width();
    vector<int> entrySizes(commSize);
    mpi::AllGather( &numLocalEntries, 1, entrySizes.data(), 1, grid.Comm() );
    vector<int> entryOffs;
    Scan( entrySizes, entryOffs );

    const auto& XDistLoc = XDist.LockedMatrix();
    if( XDistLoc.Height() == XDistLoc.LDim() )
    {
        mpi::Gather
        ( XDistLoc.LockedBuffer(), numLocalEntries,
          (T*)0, entrySizes.data(), entryOffs.data(), root, grid.Comm() );
    }
    else
    {
        vector<T> sendBuf;
        FastResize( sendBuf, numLocalEntries );
        for( Int jLoc=0; jLoc<XDistLoc.Width(); ++jLoc )
            for( Int iLoc=0; iLoc<XDistLoc.Height(); ++iLoc )
                sendBuf[iLoc+jLoc*XDistLoc.Height()] = XDistLoc(iLoc,jLoc);
        mpi::Gather
        ( sendBuf.data(), numLocalEntries,
          (T*)0, entrySizes.data(), entryOffs.data(), root, grid.Comm() );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Copy \
  ( const Matrix<T>& A, Matrix<T>& B ); \
  EL_EXTERN template void Copy \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  EL_EXTERN template void CopyFromRoot \
  ( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC>& B, bool includingViewers ); \
  EL_EXTERN template void CopyFromNonRoot \
  ( DistMatrix<T,CIRC,CIRC>& B, bool includingViewers ); \
  EL_EXTERN template void CopyFromRoot \
  ( const Matrix<T>& A, DistMatrix<T,CIRC,CIRC,BLOCK>& B, \
    bool includingViewers ); \
  EL_EXTERN template void CopyFromNonRoot \
  ( DistMatrix<T,CIRC,CIRC,BLOCK>& B, bool includingViewers ); \
  EL_EXTERN template void Copy \
  ( const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  EL_EXTERN template void Copy \
  ( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B ); \
  EL_EXTERN template void CopyFromRoot \
  ( const DistSparseMatrix<T>& ADist, SparseMatrix<T>& A ); \
  EL_EXTERN template void CopyFromNonRoot \
  ( const DistSparseMatrix<T>& ADist, int root ); \
  EL_EXTERN template void Copy \
  ( const DistMultiVec<T>& A, DistMultiVec<T>& B ); \
  EL_EXTERN template void Copy \
  ( const DistMultiVec<T>& A, AbstractDistMatrix<T>& B ); \
  EL_EXTERN template void Copy \
  ( const AbstractDistMatrix<T>& A, DistMultiVec<T>& B ); \
  EL_EXTERN template void CopyFromRoot \
  ( const DistMultiVec<T>& XDist, Matrix<T>& X ); \
  EL_EXTERN template void CopyFromNonRoot \
  ( const DistMultiVec<T>& XDist, int root );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_COPY_HPP
