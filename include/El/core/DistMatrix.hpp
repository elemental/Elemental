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

#include "./DistMatrix/Abstract.hpp"

#include "./DistMatrix/Element.hpp"
#include "./DistMatrix/Element/CIRC_CIRC.hpp"
#include "./DistMatrix/Element/MC_MR.hpp"
#include "./DistMatrix/Element/MC_STAR.hpp"
#include "./DistMatrix/Element/MD_STAR.hpp"
#include "./DistMatrix/Element/MR_MC.hpp"
#include "./DistMatrix/Element/MR_STAR.hpp"
#include "./DistMatrix/Element/STAR_MC.hpp"
#include "./DistMatrix/Element/STAR_MD.hpp"
#include "./DistMatrix/Element/STAR_MR.hpp"
#include "./DistMatrix/Element/STAR_STAR.hpp"
#include "./DistMatrix/Element/STAR_VC.hpp"
#include "./DistMatrix/Element/STAR_VR.hpp"
#include "./DistMatrix/Element/VC_STAR.hpp"
#include "./DistMatrix/Element/VR_STAR.hpp"

#include "./DistMatrix/Block.hpp"
#include "./DistMatrix/Block/CIRC_CIRC.hpp"
#include "./DistMatrix/Block/MC_MR.hpp"
#include "./DistMatrix/Block/MC_STAR.hpp"
#include "./DistMatrix/Block/MD_STAR.hpp"
#include "./DistMatrix/Block/MR_MC.hpp"
#include "./DistMatrix/Block/MR_STAR.hpp"
#include "./DistMatrix/Block/STAR_MC.hpp"
#include "./DistMatrix/Block/STAR_MD.hpp"
#include "./DistMatrix/Block/STAR_MR.hpp"
#include "./DistMatrix/Block/STAR_STAR.hpp"
#include "./DistMatrix/Block/STAR_VC.hpp"
#include "./DistMatrix/Block/STAR_VR.hpp"
#include "./DistMatrix/Block/VC_STAR.hpp"
#include "./DistMatrix/Block/VR_STAR.hpp"

namespace El {

#ifdef EL_HAVE_SCALAPACK

namespace blacs { 

template<typename scalarType>
inline int Handle( const AbstractDistMatrix<scalarType>& A )
{ return Handle( A.DistComm().comm ); }

template<typename scalarType>
inline int GridInit( int bHandle, const AbstractDistMatrix<scalarType>& A )
{
    if( A.ColDist() != MC || A.RowDist() != MR )
        LogicError
        ("Only (MC,MR) distributions are currently supported with ScaLAPACK");
    const int context =
      GridInit
      ( bHandle, A.Grid().Order()==COLUMN_MAJOR, A.ColStride(), A.RowStride() );
    DEBUG_ONLY(
      if( A.ColStride() != GridHeight(context) )
          LogicError("Grid height did not match BLACS");
      if( A.RowStride() != GridWidth(context) )
          LogicError("Grid width did not match BLACS");
      if( A.ColRank() != GridRow(context) )
          LogicError("Grid row did not match BLACS");
      if( A.RowRank() != GridCol(context) )
          LogicError("Grid col did not match BLACS");
    )
    return context;
}

} // namespace blacs

template<typename scalarType>
inline typename blacs::Desc
FillDesc( const AbstractDistMatrix<scalarType>& A, int context )
{
    if( A.ColCut() != 0 || A.RowCut() != 0 )
        LogicError("Cannot produce a meaningful descriptor if nonzero cut");
    if( A.ColDist() != MC || A.RowDist() != MR )
        LogicError
        ("Only (MC,MR) distributions are currently supported with ScaLAPACK");
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
