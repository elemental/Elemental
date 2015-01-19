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

struct DistData
{
    Dist colDist, rowDist;
    int colAlign, rowAlign; 
    int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    DistData() { }

    template<typename T>
    DistData( const AbstractDistMatrix<T>& A )
    : colDist(A.ColDist()), rowDist(A.RowDist()),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      root(A.Root()), grid(&A.Grid())
    { }
};
inline bool operator==( const DistData& A, const DistData& B )
{ return A.colDist  == B.colDist &&
         A.rowDist  == B.rowDist &&
         A.colAlign == B.colAlign &&
         A.rowAlign == B.rowAlign &&
         A.root     == B.root &&
         A.grid     == B.grid; }

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
#include "./DistMatrix/CIRC_CIRC.hpp"
#include "./DistMatrix/MC_MR.hpp"
#include "./DistMatrix/MC_STAR.hpp"
#include "./DistMatrix/MD_STAR.hpp"
#include "./DistMatrix/MR_MC.hpp"
#include "./DistMatrix/MR_STAR.hpp"
#include "./DistMatrix/STAR_MC.hpp"
#include "./DistMatrix/STAR_MD.hpp"
#include "./DistMatrix/STAR_MR.hpp"
#include "./DistMatrix/STAR_STAR.hpp"
#include "./DistMatrix/STAR_VC.hpp"
#include "./DistMatrix/STAR_VR.hpp"
#include "./DistMatrix/VC_STAR.hpp"
#include "./DistMatrix/VR_STAR.hpp"

namespace El {

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
