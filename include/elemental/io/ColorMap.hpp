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

namespace elem {

#ifdef HAVE_QT5
inline QRgb
SampleColorMap( double value, double minVal, double maxVal )
{
#ifndef RELEASE
    CallStackEntry entry("SampleColorMap");
#endif
    const double percent = (value-minVal) / (maxVal-minVal);
    const ColorMap colorMap = GetColorMap();

    int red, green, blue, alpha;
    switch( colorMap )
    {
    case RED_BLACK_GREEN:
        red = ( percent<=0.5 ? 255*(1.-2*percent) : 0 );
        green = ( percent>=0.5 ? 255*(2*(percent-0.5)) : 0 );
        blue = 0;
        alpha = 255;
        break;
    case BLUE_RED:
        red = 255*percent;
        green = 0;
        blue = 255*(1.-percent/2);
        alpha = 255;
        break;
    case GRAYSCALE:
    default:
        red = 255*percent;
        green = 255*percent;
        blue = 255*percent;
        alpha = 255;
        break;
    }

    return qRgba( red, green, blue, alpha );
}
#endif // ifdef HAVE_QT5

} // namespace elem

#endif // ifndef ELEM_IO_COLORMAP_HPP
