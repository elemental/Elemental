/*
   Copyright (c) 2009-2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
SparseLDLFactorization<Field>::SparseLDLFactorization()
{ }

template<typename Field>
void SparseLDLFactorization<Field>::Initialize
( const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE

    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );

    ldl::NestedDissection
    ( A.LockedGraph(), map_, *separator_, *info_, bisectCtrl );
    InvertMap( map_, inverseMap_ );

    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );
}

template<typename Field>
void SparseLDLFactorization<Field>::Initialize2DGridGraph
( Int gridDim0,
  Int gridDim1,
  const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE

    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );

    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, 1, A.LockedGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );

    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );
}

template<typename Field>
void SparseLDLFactorization<Field>::Initialize3DGridGraph
( Int gridDim0,
  Int gridDim1,
  Int gridDim2,
  const SparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE

    info_.reset( new ldl::NodeInfo );
    separator_.reset( new ldl::Separator );

    ldl::NaturalNestedDissection
    ( gridDim0, gridDim1, gridDim2, A.LockedGraph(),
      map_, *separator_, *info_, bisectCtrl.cutoff );
    InvertMap( map_, inverseMap_ );

    front_.reset( new ldl::Front<Field>(A,map_,*info_,hermitian) );
}

template<typename Field>
void SparseLDLFactorization<Field>::Factor( LDLFrontType frontType )
{
    EL_DEBUG_CSE
    LDL( *info_, *front_, frontType );
}

template<typename Field>
void SparseLDLFactorization<Field>::ChangeNonzeroValues
( const SparseMatrix<Field>& ANew )
{
    EL_DEBUG_CSE
    front_->Pull( ANew, map_, *info_ );
}

template<typename Field>
void SparseLDLFactorization<Field>::Solve( Matrix<Field>& B ) const
{
    EL_DEBUG_CSE
    ldl::SolveAfter( inverseMap_, *info_, *front_, B );
}

template<typename Field>
ldl::Front<Field>& SparseLDLFactorization<Field>::Front()
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
const ldl::Front<Field>& SparseLDLFactorization<Field>::Front() const
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
ldl::NodeInfo& SparseLDLFactorization<Field>::NodeInfo()
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
const ldl::NodeInfo& SparseLDLFactorization<Field>::NodeInfo() const
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
ldl::Separator& SparseLDLFactorization<Field>::Separator()
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
const ldl::Separator& SparseLDLFactorization<Field>::Separator() const
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
vector<Int>& SparseLDLFactorization<Field>::Map()
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
const vector<Int>& SparseLDLFactorization<Field>::Map() const
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
vector<Int>& SparseLDLFactorization<Field>::InverseMap()
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
const vector<Int>& SparseLDLFactorization<Field>::InverseMap() const
{
    EL_DEBUG_CSE
    return inverseMap_;
}

#define PROTO(Field) template class SparseLDLFactorization<Field>;

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
