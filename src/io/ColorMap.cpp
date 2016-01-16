/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

#ifdef EL_HAVE_QT5
QRgb SampleColorMap( double value, double minVal, double maxVal )
{
    DEBUG_ONLY(CSE cse("SampleColorMap"))
    const ColorMap colorMap = GetColorMap();
    const int numChunks = NumDiscreteColors();

    const double portion = (value-minVal) / (maxVal-minVal);
    const double discretePortion = int(portion*numChunks)/(1.*numChunks);

    int red, green, blue, alpha;
    switch( colorMap )
    {
    case RED_BLACK_GREEN:
        red = ( portion<=0.5 ? 255*(1.-2*portion) : 0 );
        green = ( portion>=0.5 ? 255*(2*(portion-0.5)) : 0 );
        blue = 0;
        alpha = 255;
        break;
    case BLUE_RED:
        red = 255*portion;
        green = 0;
        blue = 255*(1.-portion/2);
        alpha = 255;
        break;
    case GRAYSCALE_DISCRETE:
        red = 255*discretePortion;
        green = 255*discretePortion;
        blue = 255*discretePortion;
        alpha = 255;
        break;
    case GRAYSCALE:
    default:
        red = 255*portion;
        green = 255*portion;
        blue = 255*portion;
        alpha = 255;
        break;
    }

    return qRgba( red, green, blue, alpha );
}
#endif // ifdef EL_HAVE_QT5

} // namespace El
