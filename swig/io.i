/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Michael C. Grant
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

%module elem_io

%include "common.swg"
%import "elem.i"

/*
 * I/O
 */

%include "elemental/io/Print.hpp"
%include "elemental/io/Write.hpp"
%include "elemental/io/Display.hpp"
#ifdef HAVE_QT5
%include "elemental/io/Spy.hpp"
#endif

namespace elem {

OVERLOAD0_cpx(Print)
OVERLOAD1_int_UV(Print,CIRC,CIRC)
OVERLOAD1_int_UV(Print,MC,MR)
OVERLOAD1_int_UV(Print,MC,STAR)
OVERLOAD1_int_UV(Print,MD,STAR)
OVERLOAD1_int_UV(Print,MR,MC)
OVERLOAD1_int_UV(Print,MR,STAR)
OVERLOAD1_int_UV(Print,STAR,MC)
OVERLOAD1_int_UV(Print,STAR,MD)
OVERLOAD1_int_UV(Print,STAR,MR)
OVERLOAD1_int_UV(Print,STAR,STAR)
OVERLOAD1_int_UV(Print,STAR,VC)
OVERLOAD1_int_UV(Print,STAR,VR)
OVERLOAD1_int_UV(Print,VC,STAR)
OVERLOAD1_int_UV(Print,VR,STAR)

OVERLOAD0_cpx(Write)
OVERLOAD1_int_UV(Write,CIRC,CIRC)
OVERLOAD1_int_UV(Write,MC,MR)
OVERLOAD1_int_UV(Write,MC,STAR)
OVERLOAD1_int_UV(Write,MD,STAR)
OVERLOAD1_int_UV(Write,MR,MC)
OVERLOAD1_int_UV(Write,MR,STAR)
OVERLOAD1_int_UV(Write,STAR,MC)
OVERLOAD1_int_UV(Write,STAR,MD)
OVERLOAD1_int_UV(Write,STAR,MR)
OVERLOAD1_int_UV(Write,STAR,STAR)
OVERLOAD1_int_UV(Write,STAR,VC)
OVERLOAD1_int_UV(Write,STAR,VR)
OVERLOAD1_int_UV(Write,VC,STAR)
OVERLOAD1_int_UV(Write,VR,STAR)

OVERLOAD0_cpx(Display)
OVERLOAD1_int_UV(Display,CIRC,CIRC)
OVERLOAD1_int_UV(Display,MC,MR)
OVERLOAD1_int_UV(Display,MC,STAR)
OVERLOAD1_int_UV(Display,MD,STAR)
OVERLOAD1_int_UV(Display,MR,MC)
OVERLOAD1_int_UV(Display,MR,STAR)
OVERLOAD1_int_UV(Display,STAR,MC)
OVERLOAD1_int_UV(Display,STAR,MD)
OVERLOAD1_int_UV(Display,STAR,MR)
OVERLOAD1_int_UV(Display,STAR,STAR)
OVERLOAD1_int_UV(Display,STAR,VC)
OVERLOAD1_int_UV(Display,STAR,VR)
OVERLOAD1_int_UV(Display,VC,STAR)
OVERLOAD1_int_UV(Display,VR,STAR)

#ifdef HAVE_QT5
OVERLOAD0_cpx(Spy)
OVERLOAD1_int_UV(Spy,CIRC,CIRC)
OVERLOAD1_int_UV(Spy,MC,MR)
OVERLOAD1_int_UV(Spy,MC,STAR)
OVERLOAD1_int_UV(Spy,MD,STAR)
OVERLOAD1_int_UV(Spy,MR,MC)
OVERLOAD1_int_UV(Spy,MR,STAR)
OVERLOAD1_int_UV(Spy,STAR,MC)
OVERLOAD1_int_UV(Spy,STAR,MD)
OVERLOAD1_int_UV(Spy,STAR,MR)
OVERLOAD1_int_UV(Spy,STAR,STAR)
OVERLOAD1_int_UV(Spy,STAR,VC)
OVERLOAD1_int_UV(Spy,STAR,VR)
OVERLOAD1_int_UV(Spy,VC,STAR)
OVERLOAD1_int_UV(Spy,VR,STAR)
#endif

};
