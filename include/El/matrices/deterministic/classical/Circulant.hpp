/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRICES_DETERMINISTIC_CLASSICAL_CIRCULANT_HPP
#define EL_MATRICES_DETERMINISTIC_CLASSICAL_CIRCULANT_HPP

namespace El {

template<typename T> 
void Circulant( Matrix<T>& A, const Matrix<T>& a )
{
    DEBUG_ONLY(CSE cse("Circulant"))
    const Int n = a.Height();
    A.Resize( n, n );
    // NOTE: gcc (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010 segfaults here
    //       if the return type of the lambda is not manually specified.
    auto circFill = [&]( Int i, Int j ) -> T { return a.Get(Mod(i-j,n),0); };
    IndexDependentFill( A, function<T(Int,Int)>(circFill) );
}

template<typename T> 
void Circulant( Matrix<T>& A, const vector<T>& a )
{
    DEBUG_ONLY(CSE cse("Circulant"))
    const Int n = a.size();
    A.Resize( n, n );
    // NOTE: gcc (Ubuntu 5.2.1-22ubuntu2) 5.2.1 20151010 segfaults here
    //       if the return type of the lambda is not manually specified.
    auto circFill = [&]( Int i, Int j ) -> T { return a[Mod(i-j,n)]; };
    IndexDependentFill( A, function<T(Int,Int)>(circFill) );
}

template<typename T>
void Circulant( AbstractDistMatrix<T>& A, const Matrix<T>& a )
{
    DEBUG_ONLY(CSE cse("Circulant"))
    const Int n = a.Height();
    A.Resize( n, n );
    auto circFill = [&]( Int i, Int j ) -> T { return a.Get(Mod(i-j,n),0); };
    IndexDependentFill( A, function<T(Int,Int)>(circFill) );
}

template<typename T>
void Circulant( AbstractDistMatrix<T>& A, const vector<T>& a )
{
    DEBUG_ONLY(CSE cse("Circulant"))
    const Int n = a.size();
    A.Resize( n, n );
    auto circFill = [&]( Int i, Int j ) -> T { return a[Mod(i-j,n)]; };
    IndexDependentFill( A, function<T(Int,Int)>(circFill) );
}


} // namespace El

#endif // ifndef EL_MATRICES_DETERMINISTIC_CLASSICAL_CIRCULANT_HPP
