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

#include "Elemental/Matrix.hpp"
#include "Elemental/wrappers/MPI.hpp"

namespace Elemental {

// We will partially specialize for each valid distribution 
template<typename T, Distribution ColDist, Distribution RowDist> 
class DistMatrix;

} // Elemental

#include "Elemental/DistMatrix/MC_MR.hpp"
#include "Elemental/DistMatrix/MC_Star.hpp"
#include "Elemental/DistMatrix/MD_Star.hpp"
#include "Elemental/DistMatrix/MR_MC.hpp"
#include "Elemental/DistMatrix/MR_Star.hpp"
#include "Elemental/DistMatrix/Star_MC.hpp"
#include "Elemental/DistMatrix/Star_MD.hpp"
#include "Elemental/DistMatrix/Star_MR.hpp"
#include "Elemental/DistMatrix/Star_Star.hpp"
#include "Elemental/DistMatrix/Star_VC.hpp"
#include "Elemental/DistMatrix/Star_VR.hpp"
#include "Elemental/DistMatrix/VC_Star.hpp"
#include "Elemental/DistMatrix/VR_Star.hpp"

#endif /* ELEMENTAL_DISTMATRIX_HPP */
