/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_DIST_MATRIX_HPP
#define ELEMENTAL_DIST_MATRIX_HPP 1

#include "elemental/matrix.hpp"

namespace elemental {

template<typename T>
class AbstractDistMatrixBase;

template<typename T>
class AbstractDistMatrix;

template<typename T, Distribution ColDist, Distribution RowDist>
class DistMatrixBase;

template<typename T, Distribution ColDist, Distribution RowDist> 
class DistMatrix;

}

#include "elemental/dist_matrix/abstract.hpp"
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

#endif /* ELEMENTAL_DIST_MATRIX_HPP */

