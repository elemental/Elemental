/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./HermitianTridiag/L.hpp"
#include "./HermitianTridiag/LSquare.hpp"
#include "./HermitianTridiag/U.hpp"
#include "./HermitianTridiag/USquare.hpp"

#include "./HermitianTridiag/ApplyQ.hpp"

namespace El {

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiag"))
    if( uplo == LOWER )
        herm_tridiag::L( A, t );
    else
        herm_tridiag::U( A, t );
}

template<typename F> 
void HermitianTridiag
( UpperOrLower uplo, AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& tPre,
  const HermitianTridiagCtrl<F>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiag"))

    auto APtr = ReadWriteProxy<F,MC,MR>( &APre ); auto& A = *APtr;
    auto tPtr = WriteProxy<F,STAR,STAR>( &tPre ); auto& t = *tPtr;

    const Grid& g = A.Grid();
    if( ctrl.approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            herm_tridiag::L( A, t, ctrl.symvCtrl );
        else
            herm_tridiag::U( A, t, ctrl.symvCtrl );
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
        DistMatrix<F,STAR,STAR> tSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( ASquare.Participating() )
        {
            if( uplo == LOWER )
                herm_tridiag::LSquare( ASquare, tSquare, ctrl.symvCtrl );
            else
                herm_tridiag::USquare( ASquare, tSquare, ctrl.symvCtrl ); 
        }
        tSquare.MakeConsistent( true );
        A = ASquare;
        t = tSquare;

        mpi::Free( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( uplo == LOWER )
                herm_tridiag::LSquare( A, t, ctrl.symvCtrl );
            else
                herm_tridiag::USquare( A, t, ctrl.symvCtrl ); 
        }
        else
        {
            if( uplo == LOWER )
                herm_tridiag::L( A, t, ctrl.symvCtrl );
            else
                herm_tridiag::U( A, t, ctrl.symvCtrl );
        }
    }
}

namespace herm_tridiag {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("herm_tridiag::ExplicitCondensed"))
    Matrix<F> t;
    HermitianTridiag( uplo, A, t );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}


template<typename F>
void ExplicitCondensed
( UpperOrLower uplo, AbstractDistMatrix<F>& A, 
  const HermitianTridiagCtrl<F>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_tridiag::ExplicitCondensed"))
    DistMatrix<F,STAR,STAR> t(A.Grid());
    HermitianTridiag( uplo, A, t, ctrl );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

} // namespace herm_tridiag

#define PROTO(F) \
  template void HermitianTridiag \
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t ); \
  template void HermitianTridiag \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void herm_tridiag::ExplicitCondensed \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, \
    const HermitianTridiagCtrl<F>& ctrl ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, \
          AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
