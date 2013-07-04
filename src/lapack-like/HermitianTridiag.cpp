/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

#include "./HermitianTridiag/PanelL.hpp"
#include "./HermitianTridiag/PanelLSquare.hpp"
#include "./HermitianTridiag/PanelU.hpp"
#include "./HermitianTridiag/PanelUSquare.hpp"
#include "./HermitianTridiag/L.hpp"
#include "./HermitianTridiag/LSquare.hpp"
#include "./HermitianTridiag/U.hpp"
#include "./HermitianTridiag/USquare.hpp"

namespace elem {

template<typename F>
void HermitianTridiag
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianTridiag");
#endif
    if( uplo == LOWER )
        hermitian_tridiag::L( A, t );
    else
        hermitian_tridiag::U( A, t );
}

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianTridiag");
#endif
    Matrix<F> t;
    HermitianTridiag( uplo, A, t );
}

template<typename F> 
void
HermitianTridiag
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianTridiag");
#endif
    const Grid& g = A.Grid();
    const HermitianTridiagApproach approach = GetHermitianTridiagApproach();
    const GridOrder order = GetHermitianTridiagGridOrder();
    if( approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            hermitian_tridiag::L( A, t );
        else
            hermitian_tridiag::U( A, t );
    }
    else if( approach == HERMITIAN_TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        const int p = g.Size();
        const int pSqrt = int(sqrt(double(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
        if( order == COLUMN_MAJOR )
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = i+j*pSqrt;
        }
        else
        {
            for( int j=0; j<pSqrt; ++j )
                for( int i=0; i<pSqrt; ++i )
                    squareRanks[i+j*pSqrt] = j+i*pSqrt;
        }

        mpi::Group owningGroup = g.OwningGroup();
        mpi::Group squareGroup;
        mpi::GroupIncl
        ( owningGroup, squareRanks.size(), &squareRanks[0], squareGroup );

        mpi::Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt );
        DistMatrix<F> ASquare(squareGrid);
        DistMatrix<F,STAR,STAR> tSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( ASquare.Participating() )
        {
            if( uplo == LOWER )
                hermitian_tridiag::LSquare( ASquare, tSquare );
            else
                hermitian_tridiag::USquare( ASquare, tSquare ); 
        }
        tSquare.MakeConsistent();
        A = ASquare;
        t = tSquare;

        mpi::GroupFree( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( uplo == LOWER )
                hermitian_tridiag::LSquare( A, t );
            else
                hermitian_tridiag::USquare( A, t ); 
        }
        else
        {
            if( uplo == LOWER )
                hermitian_tridiag::L( A, t );
            else
                hermitian_tridiag::U( A, t );
        }
    }
}

template<typename F>
void
HermitianTridiag( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianTridiag");
#endif
    DistMatrix<F,STAR,STAR> t(A.Grid());
    HermitianTridiag( uplo, A, t );
}

#define PROTO(T) \
  template void HermitianTridiag<T>( UpperOrLower uplo, Matrix<T>& A ); \
  template void HermitianTridiag<T>( UpperOrLower uplo, DistMatrix<T>& A )

#ifndef DISABLE_FLOAT
PROTO(float);
#ifndef DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef DISABLE_FLOAT

PROTO(double);
#ifndef DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
