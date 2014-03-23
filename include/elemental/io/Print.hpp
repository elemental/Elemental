/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PRINT_HPP
#define ELEM_PRINT_HPP

namespace elem {

template<typename T>
inline void
Print( const Matrix<T>& A, std::string title="", std::ostream& os=std::cout )
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

template<typename T,Dist U,Dist V>
inline void
Print
( const DistMatrix<T,U,V>& A, std::string title="", std::ostream& os=std::cout )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( U == A.UGath && V == A.VGath )
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

template<typename T,Dist U,Dist V>
inline void
Print
( const BlockDistMatrix<T,U,V>& A, 
  std::string title="", std::ostream& os=std::cout )
{
    DEBUG_ONLY(CallStackEntry cse("Print"))
    if( U == A.UGath && V == A.VGath )
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

} // namespace elem

#endif // ifndef ELEM_PRINT_HPP
