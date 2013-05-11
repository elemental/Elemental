/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_COLORMAP_HPP
#define GRAPHICS_COLORMAP_HPP

#ifdef HAVE_QT5

namespace elem {

template<typename T>
QRgb ColorMap( T value, T minVal, T maxVal )
{
#ifndef RELEASE
    CallStackEntry entry("ColorMap");
#endif
    const double percent = double(value-minVal) / double(maxVal-minVal);

    // For now, use grey-scale
    const int red = 255*percent;
    const int green = 255*percent;
    const int blue = 255*percent;
    const int alpha = 255;

    // Red and blue mixture
    /*
    const int red = 255*percent;
    const int green = 0;
    const int blue = 255*(R(1)-percent/2);
    const int alpha = 255;
    */

    return qRgba( red, green, blue, alpha );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_COLORMAP_HPP
