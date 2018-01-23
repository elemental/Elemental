/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename T>
void ConfigurePrecision( ostream& os )
{
    // Force the full precision to be reported
    const Int numDecimals =
      BinaryToDecimalPrecision(NumMantissaBits(Base<T>()))+1;
    os.precision( numDecimals );
}

// Dense
// =====

template<typename T>
void Print( const Matrix<T>& A, string title, ostream& os )
{
    EL_DEBUG_CSE

    ostringstream msg;

    if( title != "" )
        msg << title << endl;

    ConfigurePrecision<T>( msg );

    const Int height = A.Height();
    const Int width = A.Width();
    for( Int i=0; i<height; ++i )
    {
        for( Int j=0; j<width; ++j )
            msg << A.Get(i,j) << " ";
        msg << endl;
    }
    msg << endl;
    os << msg.str();
}

template<typename T>
void Print
( const AbstractDistMatrix<T>& A, string title, ostream& os )
{
    EL_DEBUG_CSE
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Print( A.LockedMatrix(), title, os );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Print( A_CIRC_CIRC.LockedMatrix(), title, os );
    }
}

template<typename T>
void Print( const DistMultiVec<T>& X, string title, ostream& os )
{
    EL_DEBUG_CSE
    Output("Entered DistMultiVec Print with title=",title);
    if( X.Grid().Rank() == 0 )
    {
        Matrix<T> XLoc;
        CopyFromRoot( X, XLoc );
        Print( XLoc, title, os );
    }
    else
    {
        CopyFromNonRoot( X, 0 );
    }
}

void Print( const Graph& graph, string title, ostream& os )
{
    EL_DEBUG_CSE

    ostringstream msg;

    graph.AssertConsistent();
    if( title != "" )
        msg << title << endl;
    const Int numEdges = graph.NumEdges();
    const Int* srcBuf = graph.LockedSourceBuffer();
    const Int* tgtBuf = graph.LockedTargetBuffer();
    for( Int e=0; e<numEdges; ++e )
        msg << srcBuf[e] << " " << tgtBuf[e] << "\n";
    msg << endl;
    os << msg.str();
}

void Print( const DistGraph& graph, string title, ostream& os )
{
    EL_DEBUG_CSE
    graph.AssertLocallyConsistent();
    if( graph.Grid().Rank() == 0 )
    {
        Graph seqGraph;
        CopyFromRoot( graph, seqGraph );
        Print( seqGraph, title, os );
    }
    else
    {
        CopyFromNonRoot( graph, 0 );
    }
}

template<typename T>
void Print( const SparseMatrix<T>& A, string title, ostream& os )
{
    EL_DEBUG_CSE

    ostringstream msg;

    A.AssertConsistent();
    if( title != "" )
        msg << title << endl;

    ConfigurePrecision<T>( msg );

    const Int numEntries = A.NumEntries();
    const Int* srcBuf = A.LockedSourceBuffer();
    const Int* tgtBuf = A.LockedTargetBuffer();
    const T* valBuf = A.LockedValueBuffer();
    for( Int s=0; s<numEntries; ++s )
        msg << srcBuf[s] << " " << tgtBuf[s] << " " << valBuf[s] << "\n";
    msg << endl;

    os << msg.str();
}

template<typename T>
void Print( const DistSparseMatrix<T>& A, string title, ostream& os )
{
    EL_DEBUG_CSE
    A.AssertLocallyConsistent();
    if( A.Grid().Rank() == 0 )
    {
        SparseMatrix<T> ASeq;
        CopyFromRoot( A, ASeq );
        Print( ASeq, title, os );
    }
    else
    {
        CopyFromNonRoot( A, 0 );
    }
}

// Multifrontal
// ============

void PrintLocal
( const ldl::DistNodeInfo& info, string title, ostream& os )
{
    EL_DEBUG_CSE
    LogicError("This routine needs to be rewritten");
}

// Utilities
// =========

template<typename T>
void Print( const vector<T>& x, string title, ostream& os )
{
    EL_DEBUG_CSE

    ostringstream msg;

    if( title != "" )
        msg << title << endl;

    ConfigurePrecision<T>( msg );

    const Int length = x.size();
    for( Int i=0; i<length; ++i )
        msg << x[i] << " ";
    msg << endl;

    os << msg.str();
}

#define PROTO(T) \
  template void Print \
  ( const vector<T>& x, string title, ostream& os ); \
  template void Print \
  ( const Matrix<T>& A, string title, ostream& os ); \
  template void Print \
  ( const AbstractDistMatrix<T>& A, string title, ostream& os ); \
  template void Print \
  ( const DistMultiVec<T>& X, string title, ostream& os ); \
  template void Print \
  ( const SparseMatrix<T>& A, string title, ostream& os ); \
  template void Print \
  ( const DistSparseMatrix<T>& A, string title, ostream& os );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
