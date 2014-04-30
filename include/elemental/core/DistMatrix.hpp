/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DISTMATRIX_HPP
#define ELEM_DISTMATRIX_HPP

namespace elem {

struct DistData
{
    Dist colDist, rowDist;
    Int colAlign, rowAlign; 
    Int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    template<typename T,Dist U,Dist V>
    DistData( const GeneralDistMatrix<T,U,V>& A )
    : colDist(U), rowDist(V),
      colAlign(A.ColAlign()), rowAlign(A.RowAlign()),
      root(A.Root()), grid(&A.Grid())
    { }
};

} // namespace elem

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

#endif // ifndef ELEM_DISTMATRIX_HPP
