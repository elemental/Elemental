/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MAKEREAL_HPP
#define ELEM_MAKEREAL_HPP

namespace elem {

template<typename T>
inline void
MakeReal( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeReal"))
    T* ABuffer = A.Buffer();
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            ABuffer[i+j*ldim] = RealPart(ABuffer[i+j*ldim]);
}

template<typename T,Dist U,Dist V>
inline void
MakeReal( DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeReal"))
    MakeReal( A.Matrix() );
}

} // namespace elem

#endif // ifndef ELEM_MAKEREAL_HPP
