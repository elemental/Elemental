/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Trtrmm/Unblocked.hpp"
#include "./Trtrmm/LVar1.hpp"
#include "./Trtrmm/UVar1.hpp"

namespace elem {

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Trtrmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Trtrmm( Orientation orientation, UpperOrLower uplo, DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("Trtrmm");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    if( uplo == LOWER )
        internal::TrtrmmLVar1( orientation, A );
    else
        internal::TrtrmmUVar1( orientation, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
