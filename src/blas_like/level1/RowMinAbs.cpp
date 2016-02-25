/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void RowMinAbs( const Matrix<F>& A, Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbs"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    mins.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMin = limits::Max<Real>();
        for( Int j=0; j<n; ++j )
            rowMin = Min(rowMin,Abs(A.Get(i,j)));
        mins.Set( i, 0, rowMin );
    }
}

template<typename F>
void RowMinAbsNonzero
( const Matrix<F>& A, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    mins.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMin = upperBounds.Get(i,0);
        for( Int j=0; j<n; ++j )
        {
            const Real absVal = Abs(A.Get(i,j));
            if( absVal > Real(0) )
                rowMin = Min(rowMin,absVal);
        }
        mins.Set( i, 0, rowMin );
    }
}

template<typename F,Dist U,Dist V>
void RowMinAbs( const DistMatrix<F,U,V>& A, DistMatrix<Base<F>,U,STAR>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbs"))
    mins.AlignWith( A );
    mins.Resize( A.Height(), 1 );
    RowMinAbs( A.LockedMatrix(), mins.Matrix() );
    AllReduce( mins, A.RowComm(), mpi::MIN );
}

template<typename F,Dist U,Dist V>
void RowMinAbsNonzero
( const DistMatrix<F,U,V>& A, 
  const DistMatrix<Base<F>,U,STAR>& upperBounds, 
        DistMatrix<Base<F>,U,STAR>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbsNonzero"))
    if( upperBounds.ColAlign() != A.ColAlign() )
        LogicError("upperBounds was not aligned with A");
    mins.AlignWith( A );
    mins.Resize( A.Height(), 1 );
    RowMinAbsNonzero
    ( A.LockedMatrix(), upperBounds.LockedMatrix(), mins.Matrix() );
    AllReduce( mins, A.RowComm(), mpi::MIN );
}

template<typename F>
void RowMinAbs( const DistMultiVec<F>& A, DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbs"))
    mins.SetComm( A.Comm() );
    mins.Resize( A.Height(), 1 );
    RowMinAbs( A.LockedMatrix(), mins.Matrix() );
}

template<typename F>
void RowMinAbsNonzero
( const DistMultiVec<F>& A, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbsNonzero"))
    mins.SetComm( A.Comm() );
    mins.Resize( A.Height(), 1 );
    RowMinAbsNonzero
    ( A.LockedMatrix(), upperBounds.LockedMatrix(), mins.Matrix() );
}

template<typename F>
void RowMinAbs( const SparseMatrix<F>& A, Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbs"))
    typedef Base<F> Real;
    const Int m = A.Height();
    mins.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMin = limits::Max<Real>();
        const Int offset = A.RowOffset( i );
        const Int numConn = A.NumConnections( i );
        for( Int e=offset; e<offset+numConn; ++e )
            rowMin = Min(rowMin,Abs(A.Value(e)));
        mins.Set( i, 0, rowMin );
    }
}

template<typename F>
void RowMinAbsNonzero
( const SparseMatrix<F>& A, 
  const Matrix<Base<F>>& upperBounds,
        Matrix<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = A.Height();
    mins.Resize( m, 1 );
    for( Int i=0; i<m; ++i )
    {
        Real rowMin = upperBounds.Get(i,0);
        const Int offset = A.RowOffset( i );
        const Int numConn = A.NumConnections( i );
        for( Int e=offset; e<offset+numConn; ++e )
        {
            const Real absVal = Abs(A.Value(e));
            if( absVal > Real(0) )
                rowMin = Min(rowMin,absVal);
        }
        mins.Set( i, 0, rowMin );
    }
}

template<typename F>
void RowMinAbs( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbs"))
    typedef Base<F> Real;
    mins.SetComm( A.Comm() );
    mins.Resize( A.Height(), 1 );
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real rowMin = limits::Max<Real>();
        const Int offset = A.RowOffset( iLoc );
        const Int numConn = A.NumConnections( iLoc );
        for( Int e=offset; e<offset+numConn; ++e )
            rowMin = Min(rowMin,Abs(A.Value(e)));
        mins.SetLocal( iLoc, 0, rowMin ); 
    }
}

template<typename F>
void RowMinAbsNonzero
( const DistSparseMatrix<F>& A, 
  const DistMultiVec<Base<F>>& upperBounds,
        DistMultiVec<Base<F>>& mins )
{
    DEBUG_ONLY(CSE cse("RowMinAbsNonzero"))
    typedef Base<F> Real;
    mins.SetComm( A.Comm() );
    mins.Resize( A.Height(), 1 );
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        Real rowMin = upperBounds.GetLocal(iLoc,0);
        const Int offset = A.RowOffset( iLoc );
        const Int numConn = A.NumConnections( iLoc );
        for( Int e=offset; e<offset+numConn; ++e )
        {
            const Real absVal = Abs(A.Value(e));
            if( absVal > Real(0) )
                rowMin = Min(rowMin,absVal);
        }
        mins.SetLocal( iLoc, 0, rowMin ); 
    }
}

#define PROTO_DIST(F,U,V) \
  template void RowMinAbs \
  ( const DistMatrix<F,U,V>& X, DistMatrix<Base<F>,U,STAR>& mins ); \
  template void RowMinAbsNonzero \
  ( const DistMatrix<F,U,V>& X, \
    const DistMatrix<Base<F>,U,STAR>& upperBounds, \
          DistMatrix<Base<F>,U,STAR>& mins );

#define PROTO(F) \
  template void RowMinAbs \
  ( const Matrix<F>& X, Matrix<Base<F>>& mins ); \
  template void RowMinAbsNonzero \
  ( const Matrix<F>& X, \
    const Matrix<Base<F>>& upperBounds, \
          Matrix<Base<F>>& mins ); \
  template void RowMinAbs \
  ( const DistMultiVec<F>& X, DistMultiVec<Base<F>>& mins ); \
  template void RowMinAbsNonzero \
  ( const DistMultiVec<F>& X, \
    const DistMultiVec<Base<F>>& upperBounds, \
          DistMultiVec<Base<F>>& mins ); \
  template void RowMinAbs \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& mins ); \
  template void RowMinAbsNonzero \
  ( const SparseMatrix<F>& A, \
    const Matrix<Base<F>>& upperBounds, \
          Matrix<Base<F>>& mins ); \
  template void RowMinAbs \
  ( const DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& mins ); \
  template void RowMinAbsNonzero \
  ( const DistSparseMatrix<F>& A, \
    const DistMultiVec<Base<F>>& upperBounds, \
          DistMultiVec<Base<F>>& mins ); \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
