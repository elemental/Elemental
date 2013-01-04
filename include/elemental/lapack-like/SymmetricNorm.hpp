/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./SymmetricNorm/Nuclear.hpp"
#include "./SymmetricNorm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typename Base<F>::type norm = 0;
    if( type == NUCLEAR_NORM )
        norm = internal::SymmetricNuclearNorm( uplo, A );
    else if( type == TWO_NORM )
        norm = internal::SymmetricTwoNorm( uplo, A );
    else
        norm = HermitianNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typename Base<F>::type norm = 0;
    if( type == NUCLEAR_NORM )
        norm = internal::SymmetricNuclearNorm( uplo, A );
    else if( type == TWO_NORM )
        norm = internal::SymmetricTwoNorm( uplo, A );
    else
        norm = HermitianNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem
