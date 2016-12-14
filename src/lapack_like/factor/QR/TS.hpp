/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_QR_TS_HPP
#define EL_QR_TS_HPP

namespace El {
namespace qr {
namespace ts {

template<typename F>
void Reduce( const AbstractDistMatrix<F>& A, TreeData<F>& treeData )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.RowDist() != STAR )
          LogicError("Invalid row distribution for TSQR");
    )
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
    const Int logp = FlooredLog2(p);

    Matrix<F> lastZ;
    lastZ = treeData.QR0( IR(0,n), IR(0,n) );

    treeData.QRList.resize( logp );
    treeData.householderScalarsList.resize( logp );
    treeData.signatureList.resize( logp );

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
            MakeTrapezoidal( UPPER, ZTop );
            mpi::Recv( ZBot.Buffer(), n*n, partner, colComm );
        }
        else
        {
            ZBot = lastZ;
            MakeTrapezoidal( UPPER, ZBot );
            mpi::Send( ZBot.LockedBuffer(), n*n, partner, colComm );
            break;
        }

        auto& QRFact = treeData.QRList[stage];
        auto& householderScalars = treeData.householderScalarsList[stage];
        auto& signature = treeData.signatureList[stage];
        QRFact.Resize( 2*n, n, 2*n );
        householderScalars.Resize( n, 1 );
        signature.Resize( n, 1 );
        auto QRFactTop = QRFact( IR(0,n),   IR(0,n) );
        auto QRFactBot = QRFact( IR(n,2*n), IR(0,n) );
        QRFactTop = ZTop;
        QRFactBot = ZBot;

        // Note that the last QR is not performed by this routine, as many
        // higher-level routines, such as TS-SVT, are simplified if the final
        // small matrix is left alone.
        if( stage < logp-1 )
        {
            // TODO: Exploit double-triangular structure
            QR( QRFact, householderScalars, signature );
            lastZ = QRFact( IR(0,n), IR(0,n) );
        }
    }
}

template<typename F>
Matrix<F>&
RootQR( const AbstractDistMatrix<F>& A, TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root QR");
    if( p == 1 )
        return treeData.QR0;
    else
        return treeData.QRList.back();
}

template<typename F>
const Matrix<F>&
RootQR( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root QR");
    if( p == 1 )
        return treeData.QR0;
    else
        return treeData.QRList.back();
}

template<typename F>
inline Matrix<F>&
RootHouseholderScalars( const AbstractDistMatrix<F>& A, TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError
        ("This process does not have access to the root Householder scalars");
    if( p == 1 )
        return treeData.householderScalars0;
    else
        return treeData.householderScalarsList.back();
}

template<typename F>
inline const Matrix<F>&
RootHouseholderScalars
( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError
        ("This process does not have access to the root Householder scalars");
    if( p == 1 )
        return treeData.householderScalars0;
    else
        return treeData.householderScalarsList.back();
}

template<typename F>
inline Matrix<Base<F>>&
RootSignature( const AbstractDistMatrix<F>& A, TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root signature");
    if( p == 1 )
        return treeData.signature0;
    else
        return treeData.signatureList.back();
}

template<typename F>
inline const Matrix<Base<F>>&
RootSignature( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    const Int rank = mpi::Rank( A.ColComm() );
    if( rank != 0 )
        LogicError("This process does not have access to the root signature");
    if( p == 1 )
        return treeData.signature0;
    else
        return treeData.signatureList.back();
}

template<typename F>
void Scatter( AbstractDistMatrix<F>& A, const TreeData<F>& treeData )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.RowDist() != STAR )
          LogicError("Invalid row distribution for TSQR");
    )
    const Int m = A.Height();
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
    const Int logp = FlooredLog2(p);

    // Run the binary tree scatter
    Matrix<F> Z(2*n,n,2*n), ZHalf(n,n,n);
    if( rank == 0 )
        Z = RootQR( A, treeData );
    auto ZTop = Z( IR(0,n),   IR(0,n) );
    auto ZBot = Z( IR(n,2*n), IR(0,n) );
    for( Int revStage=0; revStage<logp; ++revStage )
    {
        const Int stage = (logp-1)-revStage;
        // Skip this stage if the first stage bits of our rank are not zero
        if( stage>0 && (Unsigned(rank) & ((Unsigned(1)<<stage)-1)) )
            continue;

        const Int partner = Unsigned(rank) ^ (Unsigned(1)<<stage);
        const bool top = rank < partner;
        if( top )
        {
            if( stage < logp-1 )
            {
                // Multiply by the current Q
                ZTop = ZHalf;
                Zero( ZBot );

                // TODO: Exploit sparsity?
                ApplyQ
                ( LEFT, NORMAL,
                  treeData.QRList[stage],
                  treeData.householderScalarsList[stage],
                  treeData.signatureList[stage],
                  Z );
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
    Zero( A );
    auto ATop = A.Matrix()( IR(0,n), IR(0,n) );
    ATop = ZHalf;

    // TODO: Exploit sparsity
    ApplyQ
    ( LEFT, NORMAL,
      treeData.QR0, treeData.householderScalars0, treeData.signature0,
      A.Matrix() );
}

template<typename F>
inline DistMatrix<F,STAR,STAR>
FormR( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Grid& g = A.Grid();
    DistMatrix<F,CIRC,CIRC> RRoot(g);
    if( A.ColRank() == 0 )
    {
        const Int n = A.Width();
        auto R = RootQR(A,treeData);
        auto RTop = R( IR(0,n), IR(0,n) );
        CopyFromRoot( RTop, RRoot );
        MakeTrapezoidal( UPPER, RRoot );
    }
    else
        CopyFromNonRoot( RRoot );
    DistMatrix<F,STAR,STAR> R(g);
    R = RRoot;
    return R;
}

// NOTE: This is destructive
template<typename F>
inline void
FormQ( AbstractDistMatrix<F>& A, TreeData<F>& treeData )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
    {
        A.Matrix() = treeData.QR0;
        ExpandPackedReflectors
        ( LOWER, VERTICAL, CONJUGATED, 0,
          A.Matrix(), RootHouseholderScalars(A,treeData) );
        DiagonalScale( RIGHT, NORMAL, RootSignature(A,treeData), A.Matrix() );
    }
    else
    {
        if( A.ColRank() == 0 )
        {
            ExpandPackedReflectors
            ( LOWER, VERTICAL, CONJUGATED, 0,
              RootQR(A,treeData), RootHouseholderScalars(A,treeData) );
            DiagonalScale
            ( RIGHT, NORMAL, RootSignature(A,treeData), RootQR(A,treeData) );
        }
        Scatter( A, treeData );
    }
}

} // namespace ts

template<typename F>
TreeData<F> TS( const AbstractDistMatrix<F>& A )
{
    if( A.RowDist() != STAR )
        LogicError("Invalid row distribution for TSQR");
    TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.householderScalars0, treeData.signature0 );

    const Int p = mpi::Size( A.ColComm() );
    if( p != 1 )
    {
        ts::Reduce( A, treeData );
        if( A.ColRank() == 0 )
            QR
            ( ts::RootQR(A,treeData),
              ts::RootHouseholderScalars(A,treeData),
              ts::RootSignature(A,treeData) );
    }
    return treeData;
}

template<typename F>
void ExplicitTS( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R )
{
    auto treeData = TS( A );
    Copy( ts::FormR( A, treeData ), R );
    ts::FormQ( A, treeData );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_TS_HPP
