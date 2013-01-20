/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRMM_UTIL_HPP
#define BLAS_TRMM_UTIL_HPP

namespace elem {
namespace internal {

template<typename T> 
void
SetDiagonalToOne( DistMatrix<T>& D )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int height = D.Height();
    const int localWidth = D.LocalWidth();
    const int r = D.Grid().Height();
    const int c = D.Grid().Width();
    const int colShift = D.ColShift();
    const int rowShift = D.RowShift();

    for( int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const int j = rowShift + jLoc*c;
        if( j >= 0 && j < height && (j-colShift) % r == 0 )
        {
            const int iLoc = (j-colShift)/r;
            D.SetLocal(iLoc,jLoc,1);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem


#endif // ifndef BLAS_TRMM_UTIL_HPP
