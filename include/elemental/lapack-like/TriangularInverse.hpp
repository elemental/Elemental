/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./TriangularInverse/LVar3.hpp"
#include "./TriangularInverse/UVar3.hpp"

namespace elem {

template<typename F>
inline void
TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char diagChar = UnitOrNonUnitToChar( diag );
    lapack::TriangularInverse
    ( uploChar, diagChar, A.Height(), A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
TriangularInverse
( UpperOrLower uplo, 
  UnitOrNonUnit diag, 
  DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("TriangularInverse");
#endif
    internal::TriangularInverseVar3( uplo, diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

namespace internal {

template<typename F>
inline void
TriangularInverseVar3
( UpperOrLower uplo, 
  UnitOrNonUnit diag, 
  DistMatrix<F>& A  )
{
#ifndef RELEASE
    PushCallStack("internal::TriangularInverseVar3");
#endif
    if( uplo == LOWER )
        TriangularInverseLVar3( diag, A );
    else
        TriangularInverseUVar3( diag, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal

} // namespace elem
