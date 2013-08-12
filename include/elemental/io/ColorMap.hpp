/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_COLORMAP_HPP
#define ELEM_IO_COLORMAP_HPP
#ifdef HAVE_QT5

namespace elem {

inline QRgb
ColorMap( double value, double minVal, double maxVal )
{
#ifndef RELEASE
    CallStackEntry entry("ColorMap");
#endif
    const double percent = (value-minVal) / (maxVal-minVal);

    // Grey-scale
    /*
    const Int red = 255*percent;
    const Int green = 255*percent;
    const Int blue = 255*percent;
    const Int alpha = 255;
    */

    // 0: Red, 0.5: Black, 1: Green
    const Int red = ( percent<=0.5 ? 255*(1.-2*percent) : 0 );
    const Int green = ( percent>=0.5 ? 255*(2*(percent-0.5)) : 0 );
    const Int blue = 0;
    const Int alpha = 255;

    // Red and blue mixture
    /*
    const Int red = 255*percent;
    const Int green = 0;
    const Int blue = 255*(R(1)-percent/2);
    const Int alpha = 255;
    */

    return qRgba( red, green, blue, alpha );
}

} // namespace elem

#endif // ifdef HAVE_QT5
#endif // ifndef ELEM_IO_COLORMAP_HPP
