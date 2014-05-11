/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_TS_HPP
#define ELEM_QR_TS_HPP

#include ELEM_MAKETRIANGULAR_INC
#include ELEM_EXPANDPACKEDREFLECTORS_INC
#include ELEM_QR_INC

namespace elem {
namespace qr {

template<typename F>
struct TreeData
{
    Matrix<F> QR0, t0;
    Matrix<Base<F>> d0;
    std::vector<Matrix<F>> QRList;
    std::vector<Matrix<F>> tList;
    std::vector<Matrix<Base<F>>> dList;

    TreeData( Int numStages=0 )
    : QRList(numStages), tList(numStages), dList(numStages)
    { }

    TreeData( TreeData<F>&& treeData )
    : QR0(std::move(treeData.QR0)),
      t0(std::move(treeData.t0)),
      d0(std::move(treeData.d0)),
      QRList(std::move(treeData.QRList)),
      tList(std::move(treeData.tList)),
      dList(std::move(treeData.dList))
    { }

    TreeData<F>& operator=( TreeData<F>&& treeData ) 
    {
        QR0 = std::move(treeData.QR0);
        t0 = std::move(treeData.t0);
        d0 = std::move(treeData.d0);
        QRList = std::move(treeData.QRList);
        tList = std::move(treeData.tList);
        dList = std::move(treeData.dList);
        return *this;
    }
};

namespace ts {

template<typename F,Dist U>
inline void
Reduce( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ts::Reduce"))
    const Int m =  A.Height();
    const Int n = A.Width();
    const mpi::Comm colComm = A.ColComm();
    const Int p = mpi::Size( colComm );
    if( p == 1 )
        return;
    const Int rank = mpi::Rank( colComm );
    if( m < p*n ) 
        LogicError("TSQR currently assumes height >= width*numProcesses");
    if( !PowerOfTwo(p) )
        LogicError("TSQR currently requires power-of-two number of processes");
    const Int logp = Log2(p);
    auto lastZ = LockedView( treeData.QR0, 0, 0, n, n );
    treeData.QRList.resize( logp );
    treeData.tList.resize( logp );
    treeData.dList.resize( logp );

    // Run the binary tree reduction
    Matrix<F> ZTop(n,n,n), ZBot(n,n,n);
    for( Int stage=0; stage<logp; ++stage )
    {
        // Pack, then send and receive n x n matrices
        const Int partner = Unsigned(rank) ^ (Unsigned(1)<<stage);
        const bool top = rank < partner;
        if( top )
        {
            ZTop = lastZ;
            MakeTriangular( UPPER, ZTop );
            mpi::Recv( ZBot.Buffer(), n*n, partner, colComm );
        }
        else
        {
            ZBot = lastZ;
            MakeTriangular( UPPER, ZBot );
            mpi::Send( ZBot.LockedBuffer(), n*n, partner, colComm );
            break;
        }

        auto& Q = treeData.QRList[stage];
        auto& t = treeData.tList[stage];
        auto& d = treeData.dList[stage];
        Q.Resize( 2*n, n, 2*n );
        t.Resize( n, 1 );
        d.Resize( n, 1 );
        auto QTop = View( Q, 0, 0, n, n );
        auto QBot = View( Q, n, 0, n, n );
        QTop = ZTop;
        QBot = ZBot;

        // Note that the last QR is not performed by this routine, as many
        // higher-level routines, such as TS-SVT, are simplified if the final
        // small matrix is left alone.
        if( stage < logp-1 )
        {
            // TODO: Exploit double-triangular structure
            QR( Q, t, d );
            lastZ = LockedView( Q, 0, 0, n, n );
        }
    }
}

template<typename F,Dist U>
inline Matrix<F>&
RootQR( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root QR");
    if( p == 1 )
        return treeData.QR0;
    else
        return treeData.QRList.back();
}

template<typename F,Dist U>
inline const Matrix<F>&
RootQR( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root QR");
    if( p == 1 )
        return treeData.QR0;
    else
        return treeData.QRList.back();
}

template<typename F,Dist U>
inline Matrix<F>&
RootPhases( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root phases");
    if( p == 1 )
        return treeData.t0;
    else
        return treeData.tList.back();
}

template<typename F,Dist U>
inline const Matrix<F>&
RootPhases( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root phases");
    if( p == 1 )
        return treeData.t0;
    else
        return treeData.tList.back();
}

template<typename F,Dist U>
inline Matrix<Base<F>>&
RootSignature( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root signature");
    if( p == 1 )
        return treeData.d0;
    else
        return treeData.dList.back();
}

template<typename F,Dist U>
inline const Matrix<Base<F>>&
RootSignature( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root signature");
    if( p == 1 )
        return treeData.d0;
    else
        return treeData.dList.back();
}

template<typename F,Dist U>
inline void
Scatter( DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ts::Scatter"))
    const Int m =  A.Height();
    const Int n = A.Width();
    const mpi::Comm colComm = A.ColComm();
    const Int p = mpi::Size( colComm );
    if( p == 1 )
        return;
    const Int rank = mpi::Rank( colComm );
    if( m < p*n ) 
        LogicError("TSQR currently assumes height >= width*numProcesses");
    if( !PowerOfTwo(p) )
        LogicError("TSQR currently requires power-of-two number of processes");
    const Int logp = Log2(p);

    // Run the binary tree scatter
    Matrix<F> Z(2*n,n,2*n), ZHalf(n,n,n);
    if( rank == 0 )
        Z = RootQR( A, treeData );
    auto ZTop = View( Z, 0, 0, n, n );
    auto ZBot = View( Z, n, 0, n, n );
    for( Int revStage=0; revStage<logp; ++revStage )
    {
        const Int stage = (logp-1)-revStage;
        // Skip this stage if the first stage bits of our rank are not zero
        if( stage>0 && (Unsigned(rank) & ((Unsigned(1)<<stage)-1)) )
            continue;

        const Int partner = rank ^ (1u<<stage);
        const bool top = rank < partner;
        if( top )
        {
            if( stage < logp-1 )
            {
                // Multiply by the current Q
                ZTop = ZHalf;        
                MakeZeros( ZBot );
                // TODO: Exploit sparsity?
                ApplyQ
                ( LEFT, NORMAL, 
                  treeData.QRList[stage], treeData.tList[stage], 
                  treeData.dList[stage], Z );
            }
            // Send bottom-half to partner and keep top half
            ZHalf = ZBot;
            mpi::Send( ZHalf.LockedBuffer(), n*n, partner, colComm );
            ZHalf = ZTop; 
        }
        else
        {
            // Recv top half from partner
            mpi::Recv( ZHalf.Buffer(), n*n, partner, colComm );
        }
    }

    // Apply the initial Q
    MakeZeros( A.Matrix() );
    auto ATop = View( A.Matrix(), 0, 0, n, n );
    ATop = ZHalf;
    // TODO: Exploit sparsity
    ApplyQ( LEFT, NORMAL, treeData.QR0, treeData.t0, treeData.d0, A.Matrix() );
}

template<typename F,Dist U>
inline DistMatrix<F,STAR,STAR>
FormR( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData )
{
    const Grid& g = A.Grid();
    DistMatrix<F,CIRC,CIRC> RRoot(g);
    if( A.ColRank() == 0 )
    {
        const Int n = A.Width();
        auto RTop = LockedView( RootQR(A,treeData), 0, 0, n, n );
        RRoot.CopyFromRoot( RTop );
        MakeTriangular( UPPER, RRoot );
    }
    else
        RRoot.CopyFromNonRoot();
    DistMatrix<F,STAR,STAR> R(g);
    R = RRoot;
    return R;
}

// NOTE: This is destructive
template<typename F,Dist U>
inline void
FormQ( DistMatrix<F,U,STAR>& A, TreeData<F>& treeData )
{
    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
    {
        A.Matrix() = treeData.QR0;
        ExpandPackedReflectors
        ( LOWER, VERTICAL, CONJUGATED, 0,
          A.Matrix(), RootPhases(A,treeData) );
        DiagonalScale( RIGHT, NORMAL, RootSignature(A,treeData), A.Matrix() );
    }
    else
    {
        if( A.ColRank() == 0 )
        {
            ExpandPackedReflectors
            ( LOWER, VERTICAL, CONJUGATED, 0, 
              RootQR(A,treeData), RootPhases(A,treeData) );
            DiagonalScale
            ( RIGHT, NORMAL, RootSignature(A,treeData), RootQR(A,treeData) );
        }
        Scatter( A, treeData );
    }
}

} // namespace ts

template<typename F,Dist U>
inline TreeData<F>
TS( const DistMatrix<F,U,STAR>& A )
{
    TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.t0, treeData.d0 );

    const Int p = mpi::Size( A.ColComm() );
    if( p != 1 )
    {
        ts::Reduce( A, treeData );
        if( A.ColRank() == 0 )
            QR
            ( ts::RootQR(A,treeData), ts::RootPhases(A,treeData), 
              ts::RootSignature(A,treeData) );
    }
    return treeData;
}

template<typename F,Dist U>
inline void
ExplicitTS( DistMatrix<F,U,STAR>& A, DistMatrix<F,STAR,STAR>& R )
{
    auto treeData = TS( A );
    R = ts::FormR( A, treeData );
    ts::FormQ( A, treeData );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_TS_HPP
