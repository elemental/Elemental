/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiag"))
    Matrix<F> t;
    HermitianTridiag( uplo, A, t );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

template<typename F> 
void
HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t,
  const HermitianTridiagCtrl ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiag"))
    const Grid& g = A.Grid();
    if( ctrl.approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            herm_tridiag::L( A, t );
        else
            herm_tridiag::U( A, t );
    }
    else if( ctrl.approach == HERMITIAN_TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        const Int p = g.Size();
        const Int pSqrt = Int(sqrt(double(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
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
                herm_tridiag::LSquare( ASquare, tSquare );
            else
                herm_tridiag::USquare( ASquare, tSquare ); 
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
                herm_tridiag::LSquare( A, t );
            else
                herm_tridiag::USquare( A, t ); 
        }
        else
        {
            if( uplo == LOWER )
                herm_tridiag::L( A, t );
            else
                herm_tridiag::U( A, t );
        }
    }
}

template<typename F>
void
HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A, const HermitianTridiagCtrl ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianTridiag"))
    DistMatrix<F,STAR,STAR> t(A.Grid());
    HermitianTridiag( uplo, A, t, ctrl );
    if( uplo == UPPER )
        MakeTrapezoidal( LOWER, A, 1 );
    else
        MakeTrapezoidal( UPPER, A, -1 );
}

#define PROTO(F) \
  template void HermitianTridiag\
  ( UpperOrLower uplo, Matrix<F>& A ); \
  template void HermitianTridiag\
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t ); \
  template void HermitianTridiag\
  ( UpperOrLower uplo, DistMatrix<F>& A, const HermitianTridiagCtrl ctrl ); \
  template void HermitianTridiag\
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t, \
    const HermitianTridiagCtrl ctrl ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, Matrix<F>& B ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
          DistMatrix<F>& B ); \
  template void herm_tridiag::ApplyQ \
  ( LeftOrRight side, UpperOrLower uplo, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,STAR,STAR>& t, \
          DistMatrix<F>& B );

#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
