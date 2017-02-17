/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_DISTMATRIX_HPP
#define EL_CORE_DISTMATRIX_HPP

namespace El {

template<typename DistTypeA,typename DistTypeB>
inline void AssertSameDist( const DistTypeA& distA, const DistTypeB& distB )
{
    if( distA.colDist != distB.colDist || distA.rowDist != distB.rowDist )
        RuntimeError("Matrices must have the same distribution");
}

template<typename scalarType>
class DistMultiVec;

} // namespace El

#include <El/core/DistMatrix/Abstract.hpp>

#include <El/core/DistMatrix/Element.hpp>
#include <El/core/DistMatrix/Element/CIRC_CIRC.hpp>
#include <El/core/DistMatrix/Element/MC_MR.hpp>
#include <El/core/DistMatrix/Element/MC_STAR.hpp>
#include <El/core/DistMatrix/Element/MD_STAR.hpp>
#include <El/core/DistMatrix/Element/MR_MC.hpp>
#include <El/core/DistMatrix/Element/MR_STAR.hpp>
#include <El/core/DistMatrix/Element/STAR_MC.hpp>
#include <El/core/DistMatrix/Element/STAR_MD.hpp>
#include <El/core/DistMatrix/Element/STAR_MR.hpp>
#include <El/core/DistMatrix/Element/STAR_STAR.hpp>
#include <El/core/DistMatrix/Element/STAR_VC.hpp>
#include <El/core/DistMatrix/Element/STAR_VR.hpp>
#include <El/core/DistMatrix/Element/VC_STAR.hpp>
#include <El/core/DistMatrix/Element/VR_STAR.hpp>

#include <El/core/DistMatrix/Block.hpp>
#include <El/core/DistMatrix/Block/CIRC_CIRC.hpp>
#include <El/core/DistMatrix/Block/MC_MR.hpp>
#include <El/core/DistMatrix/Block/MC_STAR.hpp>
#include <El/core/DistMatrix/Block/MD_STAR.hpp>
#include <El/core/DistMatrix/Block/MR_MC.hpp>
#include <El/core/DistMatrix/Block/MR_STAR.hpp>
#include <El/core/DistMatrix/Block/STAR_MC.hpp>
#include <El/core/DistMatrix/Block/STAR_MD.hpp>
#include <El/core/DistMatrix/Block/STAR_MR.hpp>
#include <El/core/DistMatrix/Block/STAR_STAR.hpp>
#include <El/core/DistMatrix/Block/STAR_VC.hpp>
#include <El/core/DistMatrix/Block/STAR_VR.hpp>
#include <El/core/DistMatrix/Block/VC_STAR.hpp>
#include <El/core/DistMatrix/Block/VR_STAR.hpp>

namespace El {

#ifdef EL_HAVE_SCALAPACK

namespace blacs { 

template<typename scalarType>
inline int Context( const AbstractDistMatrix<scalarType>& A )
{
    if( A.ColDist() != MC || A.RowDist() != MR )
        LogicError
        ("Only (MC,MR) distributions are currently supported with ScaLAPACK");
    return A.Grid().BlacsMCMRContext();
}

} // namespace blacs

template<typename scalarType>
inline typename blacs::Desc
FillDesc( const AbstractDistMatrix<scalarType>& A )
{
    if( A.ColCut() != 0 || A.RowCut() != 0 )
        LogicError("Cannot produce a meaningful descriptor if nonzero cut");
    const int context = blacs::Context( A );
    typename blacs::Desc desc =
        { 1, context, int(A.Height()), int(A.Width()),
          int(A.BlockHeight()), int(A.BlockWidth()),
          A.ColAlign(), A.RowAlign(), int(A.LDim()) };
    return desc;
}
#endif

template<typename scalarType>
inline void AssertSameGrids( const AbstractDistMatrix<scalarType>& A ) { }

template<typename T1,typename T2>
inline void AssertSameGrids
( const AbstractDistMatrix<T1>& A1, const AbstractDistMatrix<T2>& A2 )
{
    if( A1.Grid() != A2.Grid() )
        LogicError("Grids did not match");
}

template<typename T1,typename T2,typename... Args>
inline void AssertSameGrids
( const AbstractDistMatrix<T1>& A1, const AbstractDistMatrix<T2>& A2,
  Args&... args )
{
    if( A1.Grid() != A2.Grid() )
        LogicError("Grids did not match");
    AssertSameGrids( A2, args... );
}

template<typename scalarType>
inline void AssertSameDists( const AbstractDistMatrix<scalarType>& A ) { }

template<typename scalarType>
inline void AssertSameDists
( const AbstractDistMatrix<scalarType>& A1,
  const AbstractDistMatrix<scalarType>& A2 ) 
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
}

template<typename scalarType,typename... Args>
inline void AssertSameDists
( const AbstractDistMatrix<scalarType>& A1,
  const AbstractDistMatrix<scalarType>& A2,
  Args&... args ) 
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
    AssertSameDists( A2, args... );
}

} // namespace El

#endif // ifndef EL_CORE_DISTMATRIX_HPP
