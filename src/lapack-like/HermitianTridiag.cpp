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

template<typename R>
void HermitianTridiag( UpperOrLower uplo, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    if( uplo == LOWER )
        hermitian_tridiag::L( A );
    else
        hermitian_tridiag::U( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void HermitianTridiag
( UpperOrLower uplo, Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( uplo == LOWER )
        hermitian_tridiag::L( A, t );
    else
        hermitian_tridiag::U( A, t );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
HermitianTridiag( UpperOrLower uplo, DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    if( IsComplex<R>::val )
        throw std::logic_error("Called real routine with complex datatype");
    const Grid& g = A.Grid();
    const HermitianTridiagApproach approach = GetHermitianTridiagApproach();
    const GridOrder order = GetHermitianTridiagGridOrder();
    if( approach == HERMITIAN_TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( uplo == LOWER )
            hermitian_tridiag::L( A );
        else 
            hermitian_tridiag::U( A );
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
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<R> ASquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( uplo == LOWER )
            hermitian_tridiag::LSquare( ASquare );
        else
            hermitian_tridiag::USquare( ASquare ); 
        A = ASquare;

        mpi::GroupFree( squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( uplo == LOWER )
                hermitian_tridiag::LSquare( A );
            else
                hermitian_tridiag::USquare( A );
        }
        else
        {
            if( uplo == LOWER )
                hermitian_tridiag::L( A );
            else
                hermitian_tridiag::U( A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void
HermitianTridiag
( UpperOrLower uplo, 
  DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("HermitianTridiag");
#endif
    typedef Complex<R> C;

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
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<C> ASquare(squareGrid);
        DistMatrix<C,STAR,STAR> tSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( uplo == LOWER )
            hermitian_tridiag::LSquare( ASquare, tSquare );
        else
            hermitian_tridiag::USquare( ASquare, tSquare ); 
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
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef DISABLE_FLOAT
template void HermitianTridiag<float>( UpperOrLower uplo, Matrix<float>& A );
template void HermitianTridiag<float>( UpperOrLower uplo, DistMatrix<float>& A );
#ifndef DISABLE_COMPLEX
template void HermitianTridiag<float>( UpperOrLower uplo, Matrix<Complex<float> >& A, Matrix<Complex<float> >& t );
template void HermitianTridiag<float>( UpperOrLower uplo, DistMatrix<Complex<float> >& A, DistMatrix<Complex<float>,STAR,STAR>& t );
#endif // ifndef DISABLE_COMPLEX
#endif // ifndef DISABLE_FLOAT

template void HermitianTridiag<double>( UpperOrLower uplo, Matrix<double>& A );
template void HermitianTridiag<double>( UpperOrLower uplo, DistMatrix<double>& A );
#ifndef DISABLE_COMPLEX
template void HermitianTridiag<double>( UpperOrLower uplo, Matrix<Complex<double> >& A, Matrix<Complex<double> >& t );
template void HermitianTridiag<double>( UpperOrLower uplo, DistMatrix<Complex<double> >& A, DistMatrix<Complex<double>,STAR,STAR>& t );
#endif // ifndef DISABLE_COMPLEX

} // namespace elem
