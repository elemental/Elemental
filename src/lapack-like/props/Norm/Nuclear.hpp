/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORM_NUCLEAR_HPP
#define EL_NORM_NUCLEAR_HPP

namespace El {

template<typename F> 
Base<F> NuclearNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NuclearNorm"))
    return SchattenNorm( A, Base<F>(1) );
}

template<typename F>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianNuclearNorm"))
    return HermitianSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricNuclearNorm"))
    return SymmetricSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F,Dist U,Dist V> 
Base<F> NuclearNorm( const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("NuclearNorm"))
    return SchattenNorm( A, Base<F>(1) );
}

template<typename F,Dist U,Dist V>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianNuclearNorm"))
    return HermitianSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F,Dist U,Dist V>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricNuclearNorm"))
    return SymmetricSchattenNorm( uplo, A, Base<F>(1) );
}

} // namespace El

#endif // ifndef EL_NORM_NUCLEAR_HPP
