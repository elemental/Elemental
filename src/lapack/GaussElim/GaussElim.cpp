/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;
using namespace elemental::lapack::internal;

template<typename T>
void
elemental::lapack::GaussElim
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::GaussElim");
    if( A.GetGrid() != B.GetGrid() )
        throw logic_error( "A and B must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( A.Height() != B.Height() )
        throw logic_error( "A and B must be the same height." );
#endif
    lapack::internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        blas::Trsv( Upper, Normal, NonUnit, A, B );
    else
        blas::Trsm( Left, Upper, Normal, NonUnit, (T)1, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::GaussElim
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MC,MR>& B );

template void
elemental::lapack::GaussElim
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::GaussElim
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MC,MR>& B );

template void
elemental::lapack::GaussElim
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MC,MR>& B );
#endif

