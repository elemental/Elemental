/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_DISTMATRIX_HPP
#define EL_CORE_DISTMATRIX_HPP

namespace El {

template<typename DistTypeA,typename DistTypeB>
inline void AssertSameDist( const DistTypeA& distA, const DistTypeB& distB )
{
    if( distA.colDist != distB.colDist || distA.rowDist != distB.rowDist )
        RuntimeError("Matrices must have the same distribution");
}

template<typename T>
class DistMultiVec;

} // namespace El

#include "./DistMatrix/Abstract.hpp"

#include "./DistMatrix/Elemental.hpp"
#include "./DistMatrix/Elemental/CIRC_CIRC.hpp"
#include "./DistMatrix/Elemental/MC_MR.hpp"
#include "./DistMatrix/Elemental/MC_STAR.hpp"
#include "./DistMatrix/Elemental/MD_STAR.hpp"
#include "./DistMatrix/Elemental/MR_MC.hpp"
#include "./DistMatrix/Elemental/MR_STAR.hpp"
#include "./DistMatrix/Elemental/STAR_MC.hpp"
#include "./DistMatrix/Elemental/STAR_MD.hpp"
#include "./DistMatrix/Elemental/STAR_MR.hpp"
#include "./DistMatrix/Elemental/STAR_STAR.hpp"
#include "./DistMatrix/Elemental/STAR_VC.hpp"
#include "./DistMatrix/Elemental/STAR_VR.hpp"
#include "./DistMatrix/Elemental/VC_STAR.hpp"
#include "./DistMatrix/Elemental/VR_STAR.hpp"

#include "./DistMatrix/BlockCyclic.hpp"
#include "./DistMatrix/BlockCyclic/CIRC_CIRC.hpp"
#include "./DistMatrix/BlockCyclic/MC_MR.hpp"
#include "./DistMatrix/BlockCyclic/MC_STAR.hpp"
#include "./DistMatrix/BlockCyclic/MD_STAR.hpp"
#include "./DistMatrix/BlockCyclic/MR_MC.hpp"
#include "./DistMatrix/BlockCyclic/MR_STAR.hpp"
#include "./DistMatrix/BlockCyclic/STAR_MC.hpp"
#include "./DistMatrix/BlockCyclic/STAR_MD.hpp"
#include "./DistMatrix/BlockCyclic/STAR_MR.hpp"
#include "./DistMatrix/BlockCyclic/STAR_STAR.hpp"
#include "./DistMatrix/BlockCyclic/STAR_VC.hpp"
#include "./DistMatrix/BlockCyclic/STAR_VR.hpp"
#include "./DistMatrix/BlockCyclic/VC_STAR.hpp"
#include "./DistMatrix/BlockCyclic/VR_STAR.hpp"

namespace El {

#ifdef EL_HAVE_SCALAPACK
template<typename T>
inline typename blacs::Desc
FillDesc( const DistMatrix<T,MC,MR,BLOCK_CYCLIC>& A, int context )
{
    if( A.ColCut() != 0 || A.RowCut() != 0 )
        LogicError("Cannot produce a meaningful descriptor if nonzero cut");
    typename blacs::Desc desc =
        { 1, context, A.Height(), A.Width(), A.BlockHeight(), A.BlockWidth(),
          A.ColAlign(), A.RowAlign(), A.LDim() };
    return desc;
}
#endif

template<typename T>
inline void AssertSameGrids( const AbstractDistMatrix<T>& A ) { }

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

template<typename T>
inline void AssertSameDists( const AbstractDistMatrix<T>& A ) { }

template<typename T>
inline void AssertSameDists
( const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2 ) 
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
}

template<typename T,typename... Args>
inline void AssertSameDists
( const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2,
  Args&... args ) 
{
    if( A1.ColDist() != A2.ColDist() || A1.RowDist() != A2.RowDist() )
        LogicError("Distributions did not match");
    AssertSameDists( A2, args... );
}

} // namespace El

#endif // ifndef EL_CORE_DISTMATRIX_HPP
