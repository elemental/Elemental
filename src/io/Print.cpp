/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Print( const std::vector<T>& x, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( title != "" )
        os << title << std::endl;
    
    const Int length = x.size();
    for( Int i=0; i<length; ++i )
        os << x[i] << " ";
    os << std::endl;
}

template<typename T>
void Print( const Matrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( title != "" )
        os << title << std::endl;
    
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int i=0; i<height; ++i )
    {
        for( Int j=0; j<width; ++j )
            os << A.Get(i,j) << " ";
        os << std::endl;
    }
    os << std::endl;
}

template<typename T>
void Print
( const AbstractDistMatrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Print( A.LockedMatrix(), title, os );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Print( A_CIRC_CIRC.LockedMatrix(), title, os );
    }
}

template<typename T>
void Print
( const AbstractBlockDistMatrix<T>& A, std::string title, std::ostream& os )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Print( A.LockedMatrix(), title, os );
    }
    else
    {
        BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Print( A_CIRC_CIRC.LockedMatrix(), title, os );
    }
}

#define PROTO(T) \
  template void Print \
  ( const std::vector<T>& x, std::string title, std::ostream& os ); \
  template void Print \
  ( const Matrix<T>& A, std::string title, std::ostream& os ); \
  template void Print \
  ( const AbstractDistMatrix<T>& A, std::string title, std::ostream& os ); \
  template void Print \
  ( const AbstractBlockDistMatrix<T>& A, \
    std::string title, std::ostream& os );

#include "El/macros/Instantiate.h"

} // namespace El
