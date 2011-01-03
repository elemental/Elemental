/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

// Algorithmic controls
namespace {
lapack::internal::TridiagApproach tridiagApproach = 
    lapack::internal::TRIDIAG_NORMAL;
lapack::internal::GridOrder gridOrder = 
    lapack::internal::ROW_MAJOR;
}
void 
elemental::lapack::internal::SetTridiagApproach
( lapack::internal::TridiagApproach approach )
{ ::tridiagApproach = approach; }
void 
elemental::lapack::internal::SetTridiagSquareGridOrder
( lapack::internal::GridOrder order )
{ ::gridOrder = order; }

template<typename R>
void
elemental::lapack::Tridiag
( Shape shape, DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    const Grid& g = A.Grid();
    if( ::tridiagApproach == lapack::internal::TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( shape == Lower )
            lapack::internal::TridiagL( A );
        else 
            lapack::internal::TridiagU( A );
    }
    else if( ::tridiagApproach == lapack::internal::TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh
        int p = g.Size();
        int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
        if( ::gridOrder == lapack::internal::COL_MAJOR )
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

        MPI_Group owningGroup = g.OwningGroup();
        MPI_Group squareGroup;
        MPI_Group_incl
        ( owningGroup, squareRanks.size(), &squareRanks[0], 
          &squareGroup );

        MPI_Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<R,MC,MR> ASquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( shape == Lower )
            lapack::internal::TridiagLSquare( ASquare );
        else
            lapack::internal::TridiagUSquare( ASquare ); 
        A = ASquare;

        MPI_Group_free( &squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( shape == Lower )
                lapack::internal::TridiagLSquare( A );
            else
                lapack::internal::TridiagUSquare( A );
        }
        else
        {
            if( shape == Lower )
                lapack::internal::TridiagL( A );
            else
                lapack::internal::TridiagU( A );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<complex<R>,MC,  MR  >& A,
  DistMatrix<complex<R>,Star,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    typedef complex<R> C;

    const Grid& g = A.Grid();
    if( ::tridiagApproach == lapack::internal::TRIDIAG_NORMAL )
    {
        // Use the pipelined algorithm for nonsquare meshes
        if( shape == Lower )
            lapack::internal::TridiagL( A, t );
        else
            lapack::internal::TridiagU( A, t );
    }
    else if( ::tridiagApproach == lapack::internal::TRIDIAG_SQUARE )
    {
        // Drop down to a square mesh 
        int p = g.Size();
        int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

        std::vector<int> squareRanks(pSqrt*pSqrt);
        if( ::gridOrder == lapack::internal::COL_MAJOR )
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

        MPI_Group owningGroup = g.OwningGroup();
        MPI_Group squareGroup;
        MPI_Group_incl
        ( owningGroup, squareRanks.size(), &squareRanks[0], 
          &squareGroup );

        MPI_Comm viewingComm = g.ViewingComm();
        const Grid squareGrid( viewingComm, squareGroup, pSqrt, pSqrt );
        DistMatrix<C,MC,MR> ASquare(squareGrid);
        DistMatrix<C,Star,Star> tSquare(squareGrid);

        // Perform the fast tridiagonalization on the square grid
        ASquare = A;
        if( shape == Lower )
            lapack::internal::TridiagLSquare( ASquare, tSquare );
        else
            lapack::internal::TridiagUSquare( ASquare, tSquare ); 
        A = ASquare;
        t = tSquare;

        MPI_Group_free( &squareGroup );
    }
    else
    {
        // Use the normal approach unless we're already on a square 
        // grid, in which case we use the fast square method.
        if( g.Height() == g.Width() )
        {
            if( shape == Lower )
                lapack::internal::TridiagLSquare( A, t );
            else
                lapack::internal::TridiagUSquare( A, t ); 
        }
        else
        {
            if( shape == Lower )
                lapack::internal::TridiagL( A, t );
            else
                lapack::internal::TridiagU( A, t );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<scomplex,MC,  MR  >& A,
  DistMatrix<scomplex,Star,Star>& t );

template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<dcomplex,MC,  MR  >& A,
  DistMatrix<dcomplex,Star,Star>& t );
#endif

