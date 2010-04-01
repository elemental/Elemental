/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalLAPACKInternal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::LAPACK::Internal;

template<typename T>
void
Elemental::LAPACK::GaussElim
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("LAPACK::GaussElim");
    if( A.GetGrid() != B.GetGrid() )
        throw "A and B must be distributed over the same grid.";
    if( A.Height() != A.Width() )
        throw "A must be square.";
    if( A.Height() != B.Height() )
        throw "A and B must be the same height.";
#endif
    LAPACK::Internal::ReduceToRowEchelon( A, B );
    if( B.Width() == 1 )
        BLAS::Trsv( Upper, Normal, NonUnit, A, B );
    else
        BLAS::Trsm( Left, Upper, Normal, NonUnit, (T)1, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
Elemental::LAPACK::GaussElim
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MC,MR>& B );

template void
Elemental::LAPACK::GaussElim
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MC,MR>& B );

#ifndef WITHOUT_COMPLEX
template void
Elemental::LAPACK::GaussElim
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MC,MR>& B );

template void
Elemental::LAPACK::GaussElim
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MC,MR>& B );
#endif

