/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#include "./HermitianTridiag/HermitianPanelTridiagL.hpp"
#include "./HermitianTridiag/HermitianPanelTridiagLSquare.hpp"
#include "./HermitianTridiag/HermitianPanelTridiagU.hpp"
#include "./HermitianTridiag/HermitianPanelTridiagUSquare.hpp"
#include "./HermitianTridiag/HermitianTridiagL.hpp"
#include "./HermitianTridiag/HermitianTridiagLSquare.hpp"
#include "./HermitianTridiag/HermitianTridiagU.hpp"
#include "./HermitianTridiag/HermitianTridiagUSquare.hpp"
#include "./HermitianTridiag/LocalHermitianTridiag.hpp"

namespace elem {

template<typename R>
inline void
HermitianTridiag( UpperOrLower uplo, DistMatrix<R,MC,MR>& A )
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
        const int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

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
        DistMatrix<R,MC,MR> ASquare(squareGrid);

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
  DistMatrix<Complex<R>,MC,  MR  >& A,
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
        const int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

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
        DistMatrix<C,MC,MR> ASquare(squareGrid);
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
