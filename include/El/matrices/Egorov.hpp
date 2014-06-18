/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EGOROV_HPP
#define EL_EGOROV_HPP

namespace El {

template<typename Real,class RealFunctor>
inline void
Egorov( Matrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
             const Real theta = phase(i,j);
             return Complex<Real>(Cos(theta),Sin(theta));
         } );
}

template<typename Real,class RealFunctor>
inline void
Egorov( AbstractDistMatrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
             const Real theta = phase(i,j);
             return Complex<Real>(Cos(theta),Sin(theta));
         } );
}

template<typename Real,class RealFunctor>
inline void
Egorov
( AbstractBlockDistMatrix<Complex<Real>>& A, const RealFunctor& phase, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Egorov"))
    A.Resize( n, n );
    IndexDependentFill
    ( A, [&]( Int i, Int j )
         {
             const Real theta = phase(i,j);
             return Complex<Real>(Cos(theta),Sin(theta));
         } );
}

} // namespace El

#endif // ifndef EL_EGOROV_HPP
