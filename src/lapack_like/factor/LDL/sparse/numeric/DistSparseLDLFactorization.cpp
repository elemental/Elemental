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
DistSparseLDLFactorization<Field>::DistSparseLDLFactorization
( const DistSparseMatrix<Field>& A,
        bool hermitian,
  const BisectCtrl& bisectCtrl )
{
    EL_DEBUG_CSE

    info_.reset( new ldl::DistNodeInfo(A.Grid()) );
    separator_.reset( new ldl::DistSeparator );

    ldl::NestedDissection
    ( A.LockedDistGraph(), map_, *separator_, *info_, bisectCtrl );
    InvertMap( map_, inverseMap_ );

    front_.reset
    ( new ldl::DistFront<Field>(A,map_,*separator_,*info_,hermitian) );
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Factor( LDLFrontType frontType )
{
    EL_DEBUG_CSE
    LDL( *info_, *front_, frontType );
}

template<typename Field>
void DistSparseLDLFactorization<Field>::ChangeNonzeroValues
( const DistSparseMatrix<Field>& ANew )
{
    EL_DEBUG_CSE
    if( !formedPullMetadata_ )
    {
        ANew.MappedSources( map_, mappedSources_ );
        ANew.MappedTargets( map_, mappedTargets_, columnOffsets_ );
        formedPullMetadata_ = true;
    }
    front_->Pull
    ( ANew, map_, *separator_, *info_,
      mappedSources_, mappedTargets_, columnOffsets_ );
}

template<typename Field>
void DistSparseLDLFactorization<Field>::Solve( DistMultiVec<Field>& B ) const
{
    EL_DEBUG_CSE
    ldl::SolveAfter( inverseMap_, *info_, *front_, B );
}

template<typename Field>
ldl::DistFront<Field>& DistSparseLDLFactorization<Field>::Front()
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
const ldl::DistFront<Field>& DistSparseLDLFactorization<Field>::Front() const
{
    EL_DEBUG_CSE
    return *front_;
}

template<typename Field>
ldl::DistNodeInfo& DistSparseLDLFactorization<Field>::NodeInfo()
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
const ldl::DistNodeInfo& DistSparseLDLFactorization<Field>::NodeInfo() const
{
    EL_DEBUG_CSE
    return *info_;
}

template<typename Field>
ldl::DistSeparator& DistSparseLDLFactorization<Field>::Separator()
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
const ldl::DistSeparator& DistSparseLDLFactorization<Field>::Separator() const
{
    EL_DEBUG_CSE
    return *separator_;
}

template<typename Field>
DistMap& DistSparseLDLFactorization<Field>::Map()
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
const DistMap& DistSparseLDLFactorization<Field>::Map() const
{
    EL_DEBUG_CSE
    return map_;
}

template<typename Field>
DistMap& DistSparseLDLFactorization<Field>::InverseMap()
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
const DistMap& DistSparseLDLFactorization<Field>::InverseMap() const
{
    EL_DEBUG_CSE
    return inverseMap_;
}

template<typename Field>
ldl::DistMultiVecNodeMeta&
DistSparseLDLFactorization<Field>::DistMultiVecNodeMeta() const
{
    EL_DEBUG_CSE
    return dmvMeta_;
}

#define PROTO(Field) template class DistSparseLDLFactorization<Field>;

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
