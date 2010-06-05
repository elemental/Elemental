/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_HPP
#define ELEMENTAL_DIST_MATRIX_HPP 1

#include "elemental/matrix.hpp"
#include "elemental/wrappers/mpi.hpp"

namespace elemental {

template<typename T>
class AbstractDistMatrixBase;

template<typename T>
class AbstractDistMatrix;

template<typename T, Distribution ColDist, Distribution RowDist>
class DistMatrixBase;

template<typename T, Distribution ColDist, Distribution RowDist> 
class DistMatrix;

} // elemental

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

