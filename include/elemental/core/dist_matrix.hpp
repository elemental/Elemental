/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_HPP
#define ELEM_CORE_DISTMATRIX_HPP

namespace elem {

struct DistData
{
    Distribution colDist, rowDist;
    Int colAlign, rowAlign; 
    Int root;  // relevant for [o ,o ]/[MD,* ]/[* ,MD]
    const Grid* grid;

    template<typename T,Distribution U,Distribution V>
    DistData( const DistMatrix<T,U,V>& A )
    {
        colDist = U;
        rowDist = V;
        colAlign = A.ColAlign();
        rowAlign = A.RowAlign();
        root = A.Root();
        grid = &A.Grid();
    }
};

} // namespace elem

#include "./dist_matrix/abstract.hpp"
#include "./dist_matrix/circ_circ.hpp"
#include "./dist_matrix/mc_mr.hpp"
#include "./dist_matrix/mc_star.hpp"
#include "./dist_matrix/md_star.hpp"
#include "./dist_matrix/mr_mc.hpp"
#include "./dist_matrix/mr_star.hpp"
#include "./dist_matrix/star_mc.hpp"
#include "./dist_matrix/star_md.hpp"
#include "./dist_matrix/star_mr.hpp"
#include "./dist_matrix/star_star.hpp"
#include "./dist_matrix/star_vc.hpp"
#include "./dist_matrix/star_vr.hpp"
#include "./dist_matrix/vc_star.hpp"
#include "./dist_matrix/vr_star.hpp"

#endif // ifndef ELEM_CORE_DISTMATRIX_HPP
