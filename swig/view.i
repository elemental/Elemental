/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_view

%include "common.swg"
%import "elem.i"

%include "elemental/core/view/decl.hpp"
%include "elemental/core/partition/decl.hpp"

namespace elem {

OVERLOAD_VIEW(View)
OVERLOAD_VIEW_UV(View,CIRC,CIRC)
OVERLOAD_VIEW_UV(View,MC,STAR)
OVERLOAD_VIEW_UV(View,MR,MC)
OVERLOAD_VIEW_UV(View,MR,STAR)
OVERLOAD_VIEW_UV(View,STAR,MC)
OVERLOAD_VIEW_UV(View,STAR,MR)
OVERLOAD_VIEW_UV(View,STAR,STAR)
OVERLOAD_VIEW_UV(View,STAR,VC)
OVERLOAD_VIEW_UV(View,STAR,VR)
OVERLOAD_VIEW_UV(View,VC,STAR)
OVERLOAD_VIEW_UV(View,VR,STAR)
OVERLOAD_VIEW(LockedView)
OVERLOAD_VIEW(View1x2)
OVERLOAD_VIEW(LockedView1x2)
OVERLOAD_VIEW(View2x1)
OVERLOAD_VIEW(LockedView2x1)
OVERLOAD_VIEW(View2x2)
OVERLOAD_VIEW(LockedView2x2)

OVERLOAD_VIEW(PartitionUp)
OVERLOAD_VIEW(LockedPartitionUp)
OVERLOAD_VIEW(PartitionDown)
OVERLOAD_VIEW(LockedPartitionDown)
OVERLOAD_VIEW(PartitionLeft)
OVERLOAD_VIEW(LockedPartitionLeft)
OVERLOAD_VIEW(PartitionRight)
OVERLOAD_VIEW(LockedPartitionRight)
OVERLOAD_VIEW(PartitionUpDiagonal)
OVERLOAD_VIEW(LockedPartitionUpDiagonal)
OVERLOAD_VIEW(PartitionUpOffsetDiagonal)
OVERLOAD_VIEW(LockedPartitionUpOffsetDiagonal)
OVERLOAD_VIEW(PartitionDownDiagonal)
OVERLOAD_VIEW(LockedPartitionDownDiagonal)
OVERLOAD_VIEW(PartitionDownOffsetDiagonal)
OVERLOAD_VIEW(LockedPartitionDownOffsetDiagonal)

};
