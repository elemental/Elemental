/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANTRIDIAG_HPP
#define LAPACK_HERMITIANTRIDIAG_HPP

#include "elemental/lapack-like/Reflector.hpp"

#include "./HermitianTridiag/PanelL.hpp"
#include "./HermitianTridiag/PanelLSquare.hpp"
#include "./HermitianTridiag/PanelU.hpp"
#include "./HermitianTridiag/PanelUSquare.hpp"
#include "./HermitianTridiag/L.hpp"
#include "./HermitianTridiag/LSquare.hpp"
#include "./HermitianTridiag/U.hpp"
#include "./HermitianTridiag/USquare.hpp"
#include "./HermitianTridiag/Local.hpp"

namespace elem {

template<typename R>
inline void
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
            internal::HermitianTridiagL( A );
        else 
            internal::HermitianTridiagU( A );
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
            internal::HermitianTridiagLSquare( ASquare );
        else
            internal::HermitianTridiagUSquare( ASquare ); 
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
                internal::HermitianTridiagLSquare( A );
            else
                internal::HermitianTridiagUSquare( A );
        }
        else
        {
            if( uplo == LOWER )
                internal::HermitianTridiagL( A );
            else
                internal::HermitianTridiagU( A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
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
            internal::HermitianTridiagL( A, t );
        else
            internal::HermitianTridiagU( A, t );
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
            internal::HermitianTridiagLSquare( ASquare, tSquare );
        else
            internal::HermitianTridiagUSquare( ASquare, tSquare ); 
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
                internal::HermitianTridiagLSquare( A, t );
            else
                internal::HermitianTridiagUSquare( A, t ); 
        }
        else
        {
            if( uplo == LOWER )
                internal::HermitianTridiagL( A, t );
            else
                internal::HermitianTridiagU( A, t );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_HPP
