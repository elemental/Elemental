/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_DISTMATRIX_HPP
#define ELEMENTAL_DISTMATRIX_HPP 1

#include "elemental/matrix.hpp"
#include "elemental/wrappers/mpi.hpp"

namespace elemental {

// We will partially specialize for each valid distribution 
template<typename T, Distribution ColDist, Distribution RowDist> 
class DistMatrix;

} // elemental

#include "elemental/dist_matrix/mc_mr.hpp"
#include "elemental/dist_matrix/mc_star.hpp"
#include "elemental/dist_matrix/md_star.hpp"
#include "elemental/dist_matrix/mr_mc.hpp"
#include "elemental/dist_matrix/mr_star.hpp"
#include "elemental/dist_matrix/star_mc.hpp"
#include "elemental/dist_matrix/star_md.hpp"
#include "elemental/dist_matrix/star_mr.hpp"
#include "elemental/dist_matrix/star_star.hpp"
#include "elemental/dist_matrix/star_vc.hpp"
#include "elemental/dist_matrix/star_vr.hpp"
#include "elemental/dist_matrix/vc_star.hpp"
#include "elemental/dist_matrix/vr_star.hpp"

#endif /* ELEMENTAL_DISTMATRIX_HPP */
