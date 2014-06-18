/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RIEMANN_HPP
#define EL_RIEMANN_HPP

namespace El {

template<typename T> 
inline void
Riemann( Matrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riemann"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( ((j+2)%(i+2))==0 ) { return T(i+1); }
           else                   { return T(-1);  } } );
}

template<typename T>
inline void
Riemann( AbstractDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riemann"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( ((j+2)%(i+2))==0 ) { return T(i+1); }
           else                   { return T(-1);  } } );
}

template<typename T>
inline void
Riemann( AbstractBlockDistMatrix<T>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Riemann"))
    R.Resize( n, n );
    IndexDependentFill
    ( R, []( Int i, Int j ) 
         { if( ((j+2)%(i+2))==0 ) { return T(i+1); }
           else                   { return T(-1);  } } );
}

} // namespace El

#endif // ifndef EL_RIEMANN_HPP
