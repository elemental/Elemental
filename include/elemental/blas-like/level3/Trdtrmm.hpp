/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Trdtrmm/Unblocked.hpp"
#include "./Trdtrmm/LVar1.hpp"
#include "./Trdtrmm/UVar1.hpp"

namespace elem {

template<typename F>
inline void
Trdtrmm( Orientation orientation, UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Trdtrdmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( orientation, A );
    else
        internal::TrdtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
Trdtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Trdtrmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrdtrmmLVar1( orientation, A );
    else
        internal::TrdtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
