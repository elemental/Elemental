/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_HPP
#define EL_DISTMATRIX_HPP

namespace El {

struct DistData
{
    Dist colDist, rowDist;
    Int colAlign, rowAlign; 
    Int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    DistData() { }

    template<typename T,Dist U,Dist V>
    DistData( const GeneralDistMatrix<T,U,V>& A )
    : colDist(U), rowDist(V),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      root(A.Root()), grid(&A.Grid())
    { }
};

template<typename DistTypeA,typename DistTypeB>
inline void AssertSameDist( const DistTypeA& distA, const DistTypeB& distB )
{
    if( distA.colDist != distB.colDist || distA.rowDist != distB.rowDist )
        RuntimeError("Matrices must have the same distribution");
}

} // namespace El

#include "./DistMatrix/Abstract.hpp"
#include "./DistMatrix/General.hpp"
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

#endif // ifndef EL_DISTMATRIX_HPP
