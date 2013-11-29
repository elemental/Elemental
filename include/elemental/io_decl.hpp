/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_DECL_HPP
#define ELEM_IO_DECL_HPP

namespace elem {

// TODO: Distributed file formats?
namespace file_format_wrapper {
enum FileFormat
{
    ASCII,
    MATLAB_ASCII,
    BMP,
    JPG,
    JPEG,
    PNG,
    PPM,
    XBM,
    XPM
};
}
using namespace file_format_wrapper;

// TODO: Many more color maps
namespace color_map_wrapper {
enum ColorMap
{
    GRAYSCALE,
    RED_BLACK_GREEN,
    BLUE_RED
};
}
using namespace color_map_wrapper;

void SetColorMap( ColorMap colorMap );
ColorMap GetColorMap();

// When Elemental is finalized, if no window was opened, then it must call 
// app.exit() instead
void OpenedWindow();

// For keeping track of the extreme values visualized so far
double MinRealWindowVal();
double MaxRealWindowVal();
double MinImagWindowVal();
double MaxImagWindowVal();
void UpdateMinRealWindowVal( double minVal );
void UpdateMaxRealWindowVal( double maxVal );
void UpdateMinImagWindowVal( double minVal );
void UpdateMaxImagWindowVal( double maxVal );

} // namespace elem

#ifdef HAVE_QT5
# include "elemental/io/display_widget_decl.hpp"
# include "elemental/io/display_window_decl.hpp"
# include "elemental/io/complex_display_window_decl.hpp"
#endif // ifdef HAVE_QT5

#endif // ifndef ELEM_IO_DECL_HPP
