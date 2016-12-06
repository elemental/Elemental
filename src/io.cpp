/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace {
using namespace El;

ColorMap colorMap=RED_BLACK_GREEN;
Int numDiscreteColors = 15;

}

namespace El {

void SetColorMap( ColorMap map )
{ ::colorMap = map; }

ColorMap GetColorMap()
{ return ::colorMap; }

void SetNumDiscreteColors( Int numChunks )
{ ::numDiscreteColors = numChunks; }

Int NumDiscreteColors()
{ return ::numDiscreteColors; }

} // namespace El
