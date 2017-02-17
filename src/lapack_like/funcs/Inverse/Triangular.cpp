/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Triangular/LVar3.hpp"
#include "./Triangular/UVar3.hpp"

namespace El {
namespace triang_inv {

template<typename Field>
void Var3( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<Field>& A  )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
}

template<typename Field>
void Var3
( UpperOrLower uplo, UnitOrNonUnit diag, AbstractDistMatrix<Field>& A  )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        LVar3( diag, A );
    else
        UVar3( diag, A );
}

} // namespace triang_inv

template<typename Field>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<Field>& A )
{
    EL_DEBUG_CSE
    triang_inv::Var3( uplo, diag, A );
}

template<typename Field>
void TriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, AbstractDistMatrix<Field>& A  )
{
    EL_DEBUG_CSE
    triang_inv::Var3( uplo, diag, A );
}

template<typename Field>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<Field,STAR,STAR>& A )
{
    EL_DEBUG_CSE
    TriangularInverse( uplo, diag, A.Matrix() );
}

#define PROTO(Field) \
  template void TriangularInverse \
  ( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<Field>& A ); \
  template void TriangularInverse \
  ( UpperOrLower uplo, UnitOrNonUnit diag, AbstractDistMatrix<Field>& A ); \
  template void LocalTriangularInverse \
  ( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<Field,STAR,STAR>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
