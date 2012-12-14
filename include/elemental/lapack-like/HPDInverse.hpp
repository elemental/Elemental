/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./HPDInverse/LVar2.hpp"
#include "./HPDInverse/UVar2.hpp"

namespace elem {

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("HPDInverse");
#endif
    if( uplo == LOWER )
        internal::HPDInverseLVar2( A );
    else
        internal::HPDInverseUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
