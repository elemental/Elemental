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
    Int colAlignment, rowAlignment; 
    Int root;     // only relevant for [o ,o ]
    Int diagPath; // only relevant for [MD,* ]/[* ,MD] distributions
    const Grid* grid;
};
} // namespace elem

#include "elemental/core/dist_matrix/abstract.hpp"
#include "elemental/core/dist_matrix/circ_circ.hpp"
#include "elemental/core/dist_matrix/mc_mr.hpp"
#include "elemental/core/dist_matrix/mc_star.hpp"
#include "elemental/core/dist_matrix/md_star.hpp"
#include "elemental/core/dist_matrix/mr_mc.hpp"
#include "elemental/core/dist_matrix/mr_star.hpp"
#include "elemental/core/dist_matrix/star_mc.hpp"
#include "elemental/core/dist_matrix/star_md.hpp"
#include "elemental/core/dist_matrix/star_mr.hpp"
#include "elemental/core/dist_matrix/star_star.hpp"
#include "elemental/core/dist_matrix/star_vc.hpp"
#include "elemental/core/dist_matrix/star_vr.hpp"
#include "elemental/core/dist_matrix/vc_star.hpp"
#include "elemental/core/dist_matrix/vr_star.hpp"

#endif // ifndef ELEM_CORE_DISTMATRIX_HPP
