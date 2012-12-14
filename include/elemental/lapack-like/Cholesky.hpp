/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Cholesky/LVar2.hpp"
#include "./Cholesky/LVar3.hpp"
#include "./Cholesky/LVar3Square.hpp"
#include "./Cholesky/UVar2.hpp"
#include "./Cholesky/UVar3.hpp"
#include "./Cholesky/UVar3Square.hpp"

namespace elem {

template<typename F>
inline void
Cholesky( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    lapack::Cholesky( uploChar, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Cholesky( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Cholesky");
#endif
    const Grid& g = A.Grid();

    // TODO: Come up with a better routing mechanism
    if( g.Height() == g.Width() )
    {
        if( uplo == LOWER )
            internal::CholeskyLVar3Square( A );
        else
            internal::CholeskyUVar3Square( A );
    }
    else
    {
        if( uplo == LOWER )
            internal::CholeskyLVar3( A );
        else
            internal::CholeskyUVar3( A );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
