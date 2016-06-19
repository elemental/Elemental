/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace herm_tridiag {

template<typename F>
void Ger2Sub
( const DistMatrix<F,MC,STAR>& x_MC,
  const DistMatrix<F,MC,STAR>& y_MC,
  const DistMatrix<F,MR,STAR>& x_MR,
  const DistMatrix<F,MR,STAR>& y_MR,
        DistMatrix<F,MC,MR  >& A )
{
    DEBUG_CSE
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const F* x_MC_Buf = x_MC.LockedBuffer();
    const F* x_MR_Buf = x_MR.LockedBuffer();
    const F* y_MC_Buf = y_MC.LockedBuffer();
    const F* y_MR_Buf = y_MR.LockedBuffer();

    Matrix<F>& ALoc = A.Matrix();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const F delta = Conj(x_MR_Buf[jLoc]);
        const F gamma = Conj(y_MR_Buf[jLoc]);
        F* aLoc = &ALoc(0,jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            aLoc[iLoc] -= y_MC_Buf[iLoc]*delta + x_MC_Buf[iLoc]*gamma;
    }
}

} // namespace herm_tridiag
} // namespace El

#include "./HermitianTridiag/L.hpp"
#include "./HermitianTridiag/LSquare.hpp"
#include "./HermitianTridiag/U.hpp"
#include "./HermitianTridiag/USquare.hpp"

#include "./HermitianTridiag/ApplyQ.hpp"

namespace El {

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& phase )
{
    DEBUG_CSE
    if( uplo == LOWER )
        herm_tridiag::L( A, phase );
    else
        herm_tridiag::U( A, phase );
}

template<typename F> 
void HermitianTridiag
( UpperOrLower uplo,
  ElementalMatrix<F>& APre,
  ElementalMatrix<F>& phasePre,
  const HermitianTridiagCtrl<F>& ctrl )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR> phaseProx( phasePre );
    auto& A = AProx.Get();
    auto& phase = phaseProx.Get();

    const Grid& g = A.Grid();
    if( ctrl.approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            herm_tridiag::L( A, phase, ctrl.symvCtrl );
        else
            herm_tridiag::U( A, phase, ctrl.symvCtrl );
    }
    else if( ctrl.approach == HERMITIAN_TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        const Int p = g.Size();
        const Int pSqrt = Int(sqrt(double(p)));

        vector<int> squareRanks(pSqrt*pSqrt);
        if( ctrl.order == g.Order() )
        {
            for( Int j=0; j<pSqrt; ++j )
                for( Int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = i+j*pSqrt;
        }
        else
        {
            for( Int j=0; j<pSqrt; ++j )
                for( Int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = j+i*pSqrt;
        }

        mpi::Group owningGroup = g.OwningGroup();
        mpi::Group squareGroup;
        mpi::Incl
        ( owningGroup, squareRanks.size(), squareRanks.data(), squareGroup );

        mpi::Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt );
        DistMatrix<F> ASquare(squareGrid);
        DistMatrix<F,STAR,STAR> phaseSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( ASquare.Participating() )
        {
            if( uplo == LOWER )
                herm_tridiag::LSquare( ASquare, phaseSquare, ctrl.symvCtrl );
            else
                herm_tridiag::USquare( ASquare, phaseSquare, ctrl.symvCtrl );
        }
        phaseSquare.MakeConsistent( true );
        A = ASquare;
        phase = phaseSquare;

        mpi::Free( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( uplo == LOWER )
                herm_tridiag::LSquare( A, phase, ctrl.symvCtrl );
            else
                herm_tridiag::USquare( A, phase, ctrl.symvCtrl ); 
        }
        else
        {
            if( uplo == LOWER )
                herm_tridiag::L( A, phase, ctrl.symvCtrl );
            else
                herm_tridiag::U( A, phase, ctrl.symvCtrl );
        }
    }
}

namespace herm_tridiag {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> phase;
    HermitianTridiag( uplo, A, phase );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F>
void ExplicitCondensed
( UpperOrLower uplo,
  ElementalMatrix<F>& A, 
  const HermitianTridiagCtrl<F>& ctrl )
{
    DEBUG_CSE
    DistMatrix<F,STAR,STAR> phase(A.Grid());
    HermitianTridiag( uplo, A, phase, ctrl );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

} // namespace herm_tridiag

#define PROTO(F) \
  template void HermitianTridiag \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<F>& phase ); \
  template void HermitianTridiag \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<F>& phase, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& phase, \
          Matrix<F>& B ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& phase, \
          ElementalMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
