/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_QT5
# include "El/io/DisplayWindow-premoc.hpp"
# include "El/io/ComplexDisplayWindow-premoc.hpp"
# include <QApplication>
#endif

namespace El {

void ProcessEvents( int numMsecs )
{
#ifdef EL_HAVE_QT5
    QCoreApplication::instance()->processEvents
    ( QEventLoop::AllEvents, numMsecs );
#endif
}

template<typename Real>
void Display( const Matrix<Real>& A, string title )
{
    DEBUG_CSE
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    // Convert A to double-precision since Qt's MOC does not support templates
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<double>* ADouble = new Matrix<double>( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ADouble->Set( i, j, double(A.Get(i,j)) );

    QString qTitle = QString::fromStdString( title );
    DisplayWindow* displayWindow = new DisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename Real>
void Display( const Matrix<Complex<Real>>& A, string title )
{
    DEBUG_CSE
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
    {
        Print( A, title );
        return;
    }

    // Convert A to double-precision since Qt's MOC does not support templates
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Complex<double>>* ADouble = new Matrix<Complex<double>>( m, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Complex<Real> alpha = A.Get(i,j);
            const Complex<double> alphaDouble = 
                Complex<double>(alpha.real(),alpha.imag()); 
            ADouble->Set( i, j, alphaDouble );
        }
    }

    QString qTitle = QString::fromStdString( title );
    ComplexDisplayWindow* displayWindow = new ComplexDisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T>
void Display( const AbstractDistMatrix<T>& A, string title )
{
    DEBUG_CSE
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Display( A.LockedMatrix(), title );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Display( A_CIRC_CIRC.Matrix(), title );
    }
}

template<typename T>
void Display( const DistMultiVec<T>& X, string title )
{
    DEBUG_CSE
    const int commRank = mpi::Rank( X.Comm() );
    if( commRank == 0 )
    {
        Matrix<T> XLoc;
        CopyFromRoot( X, XLoc );
        Display( XLoc, title );
    }
    else
    {
        CopyFromNonRoot( X, 0 );
    }
}

void Display( const Graph& graph, string title )
{
    DEBUG_CSE
#ifdef HAVE_QT5
    graph.AssertConsistent();
    auto graphMat = new Matrix<int>;
    const int m = graph.NumTargets();
    const int n = graph.NumSources();
    Zeros( *graphMat, m, n );

    const int numEdges = graph.NumEdges();
    const int* srcBuf = graph.LockedSourceBuffer();
    const int* tgtBuf = graph.LockedTargetBuffer();
    for( int e=0; e<numEdges; ++e )
        graphMat->Set( tgtBuf[e], srcBuf[e], 1 );

    QString qTitle = QString::fromStdString( title );
    auto spyWindow = new SpyWindow;
    spyWindow->Spy( graphMat, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( graph, title );
#endif
}

void Display( const DistGraph& graph, string title )
{
    DEBUG_CSE
    const mpi::Comm comm = graph.Comm();
    const int commRank = mpi::Rank( comm );

    if( commRank == 0 )
    {
        Graph seqGraph;
        CopyFromRoot( graph, seqGraph );
        Display( seqGraph, title );
    }
    else
    {
        CopyFromNonRoot( graph, 0 );
    }
}

template<typename Real>
void Display( const SparseMatrix<Real>& A, string title )
{
    DEBUG_CSE
#ifdef HAVE_QT5
    A.AssertConsistent();
    auto AFull = new Matrix<double>;
    const int m = A.Height();
    const int n = A.Width();
    Zeros( *AFull, m, n );

    const int numEntries = A.NumEntries();
    const int* srcBuf = A.LockedSourceBuffer();
    const int* tgtBuf = A.LockedTargetBuffer();
    const Real* valBuf = A.LockedValueBuffer();
    for( int s=0; s<numEntries; ++s )
        AFull->Set( tgtBuf[s], srcBuf[s], double(valBuf[s]) );

    QString qTitle = QString::fromStdString( title );
    auto displayWindow = new DisplayWindow;
    displayWindow->Display( AFull, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename Real>
void Display( const SparseMatrix<Complex<Real>>& A, string title )
{
    DEBUG_CSE
#ifdef HAVE_QT5
    A.AssertConsistent();
    auto AFull = new Matrix<Complex<double>>;
    const int m = A.Height();
    const int n = A.Width();
    Zeros( *AFull, m, n );

    const int numEntries = A.NumEntries();
    const int* srcBuf = A.LockedSourceBuffer();
    const int* tgtBuf = A.LockedTargetBuffer();
    const Complex<Real>* valBuf = A.LockedValueBuffer();
    for( int s=0; s<numEntries; ++s )
    {
        const Complex<double> alpha =
            Complex<double>(valBuf[s].real,valBuf[s].imag);
        AFull->Set( tgtBuf[s], srcBuf[s], alpha );
    }

    QString qTitle = QString::fromStdString( title );
    auto displayWindow = new ComplexDisplayWindow;
    displayWindow->Display( AFull, qTitle );
    displayWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    Print( A, title );
#endif
}

template<typename T>
void Display( const DistSparseMatrix<T>& A, string title )
{
    DEBUG_CSE
    A.AssertLocallyConsistent();
    const mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank( comm );
    
    if( commRank == 0 )
    {
        SparseMatrix<T> ASeq; 
        CopyFromRoot( A, ASeq );
        Display( ASeq, title );
    }
    else 
    {
        CopyFromNonRoot( A, 0 );
    }
}

void DisplayLocal
( const ldl::DistNodeInfo& info, bool beforeFact, string title )
{
    DEBUG_CSE
#ifdef HAVE_QT5
    const int n = info.distNodes.back().size + info.distNodes.back().off;
    auto graphMat = new Matrix<int>;
    Zeros( *graphMat, n, n );

    const int numLocal = info.localNodes.size();
    for( int s=0; s<numLocal; ++s )
    {
        const ldl::NodeInfo& node = info.localNodes[s];
        for( int j=0; j<node.size; ++j )
            for( int i=0; i<node.size; ++i )
                graphMat->Set( i+node.off, j+node.off, 1 );
        if( beforeFact )
        {
            const int origStructSize = node.origLowerStruct.size();
            for( int i=0; i<origStructSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.origLowerStruct[i], j+node.off, 1 );
        }
        else
        {
            const int structSize = node.lowerStruct.size();
            for( int i=0; i<structSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.lowerStruct[i], j+node.off, 1 );
        }
    }

    const int numDist = info.distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        const ldl::DistNodeInfo& node = info.distNodes[s];
        for( int j=0; j<node.size; ++j )
            for( int i=0; i<node.size; ++i )
                graphMat->Set( i+node.off, j+node.off, 1 );
        if( beforeFact )
        {
            const int origStructSize = node.origLowerStruct.size();
            for( int i=0; i<origStructSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.origLowerStruct[i], j+node.off, 1 );
        }
        else
        {
            const int structSize = node.lowerStruct.size();
            for( int i=0; i<structSize; ++i )
                for( int j=0; j<node.size; ++j )
                    graphMat->Set( node.lowerStruct[i], j+node.off, 1 );
        }
    }

    QString qTitle = QString::fromStdString( title );
    auto spyWindow = new SpyWindow;
    spyWindow->Spy( graphMat, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    PrintLocal( info );
#endif
}

#define PROTO(T) \
  template void Display( const Matrix<T>& A, string title ); \
  template void Display( const AbstractDistMatrix<T>& A, string title ); \
  template void Display( const DistMultiVec<T>& X, string title ); \
  template void Display( const SparseMatrix<T>& A, string title ); \
  template void Display( const DistSparseMatrix<T>& A, string title );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
