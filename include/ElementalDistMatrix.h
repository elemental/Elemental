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
#ifndef ELEMENTAL_DISTMATRIX_H
#define ELEMENTAL_DISTMATRIX_H 1

#include "ElementalMatrix.h"
#include "wrappers/MPI.h"

namespace Elemental
{
    // We will partially specialize for each valid distribution 
    template<typename T, Distribution ColDist, Distribution RowDist> 
    class DistMatrix;
}

#include "ElementalDistMatrix_MC_MR.h"
#include "ElementalDistMatrix_MC_Star.h"
#include "ElementalDistMatrix_MD_Star.h"
#include "ElementalDistMatrix_MR_MC.h"
#include "ElementalDistMatrix_MR_Star.h"
#include "ElementalDistMatrix_Star_MC.h"
#include "ElementalDistMatrix_Star_MD.h"
#include "ElementalDistMatrix_Star_MR.h"
#include "ElementalDistMatrix_Star_Star.h"
#include "ElementalDistMatrix_Star_VC.h"
#include "ElementalDistMatrix_Star_VR.h"
#include "ElementalDistMatrix_VC_Star.h"
#include "ElementalDistMatrix_VR_Star.h"

#endif /* ELEMENTAL_DISTMATRIX_H */
