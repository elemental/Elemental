/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T,typename S>
void SetDiagonal( Matrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Set(i,j,alpha);
    }
}

template<typename T,typename S>
void SetDiagonal( AbstractDistMatrix<T>& A, S alpha, Int offset )
{
    DEBUG_ONLY(CallStackEntry cse("SetDiagonal"))
    const Int height = A.Height();
    const Int width = A.Width();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height && A.IsLocalRow(i) )
        {
            const Int iLoc = A.LocalRow(i);
            A.SetLocal( iLoc, jLoc, alpha );
        }
    }
}

#define PROTO_TYPES(T,S) \
  template void SetDiagonal( Matrix<T>& A, S alpha, Int offset ); \
  template void SetDiagonal \
  ( AbstractDistMatrix<T>& A, S alpha, Int offset ); \

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T) 

#define PROTO_CPX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T) \
  PROTO_TYPES(T,Base<T>)

PROTO_INT(Int);
PROTO_REAL(float);
PROTO_REAL(double);
PROTO_CPX(Complex<float>);
PROTO_CPX(Complex<double>);

} // namespace El
