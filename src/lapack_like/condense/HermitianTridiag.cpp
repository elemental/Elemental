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
    EL_DEBUG_CSE
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

#include "./HermitianTridiag/LowerBlocked.hpp"
#include "./HermitianTridiag/LowerBlockedSquare.hpp"
#include "./HermitianTridiag/UpperBlocked.hpp"
#include "./HermitianTridiag/UpperBlockedSquare.hpp"

#include "./HermitianTridiag/ApplyQ.hpp"

namespace El {

template<typename F>
void HermitianTridiag
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& householderScalars )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        herm_tridiag::LowerBlocked( A, householderScalars );
    else
        herm_tridiag::UpperBlocked( A, householderScalars );
}

template<typename F> 
void HermitianTridiag
( UpperOrLower uplo,
  AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalarsPre,
  const HermitianTridiagCtrl<F>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& A = AProx.Get();
    auto& householderScalars = householderScalarsProx.Get();

    const Grid& grid = A.Grid();
    if( ctrl.approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            herm_tridiag::LowerBlocked( A, householderScalars, ctrl.symvCtrl );
        else
            herm_tridiag::UpperBlocked( A, householderScalars, ctrl.symvCtrl );
    }
    else if( ctrl.approach == HERMITIAN_TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        const Int p = grid.Size();
        const Int pSqrt = Int(sqrt(double(p)));

        vector<int> squareRanks(pSqrt*pSqrt);
        if( ctrl.order == grid.Order() )
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

        mpi::Group owningGroup = grid.OwningGroup();
        mpi::Group squareGroup;
        mpi::Incl
        ( owningGroup, squareRanks.size(), squareRanks.data(), squareGroup );

        mpi::Comm viewingComm = grid.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt );
        DistMatrix<F> ASquare(squareGrid);
        DistMatrix<F,STAR,STAR> householderScalarsSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( ASquare.Participating() )
        {
            if( uplo == LOWER )
                herm_tridiag::LowerBlockedSquare
                ( ASquare, householderScalarsSquare, ctrl.symvCtrl );
            else
                herm_tridiag::UpperBlockedSquare
                ( ASquare, householderScalarsSquare, ctrl.symvCtrl );
        }
        const bool includeViewers = true;
        householderScalarsSquare.MakeConsistent( includeViewers );
        A = ASquare;
        householderScalars = householderScalarsSquare;

        mpi::Free( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( grid.Height() == grid.Width() )
        {
            if( uplo == LOWER )
                herm_tridiag::LowerBlockedSquare
                ( A, householderScalars, ctrl.symvCtrl );
            else
                herm_tridiag::UpperBlockedSquare
                ( A, householderScalars, ctrl.symvCtrl ); 
        }
        else
        {
            if( uplo == LOWER )
                herm_tridiag::LowerBlocked
                ( A, householderScalars, ctrl.symvCtrl );
            else
                herm_tridiag::UpperBlocked
                ( A, householderScalars, ctrl.symvCtrl );
        }
    }
}

namespace herm_tridiag {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A )
{
    EL_DEBUG_CSE
    Matrix<F> householderScalars;
    HermitianTridiag( uplo, A, householderScalars );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F>
void ExplicitCondensed
( UpperOrLower uplo,
  AbstractDistMatrix<F>& A, 
  const HermitianTridiagCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<F,STAR,STAR> householderScalars(A.Grid());
    HermitianTridiag( uplo, A, householderScalars, ctrl );
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
    Matrix<F>& householderScalars ); \
  template void HermitianTridiag \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalars, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<F>& A, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
          Matrix<F>& B ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
          AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
