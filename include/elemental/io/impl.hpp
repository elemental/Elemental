/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_IMPL_HPP
#define ELEM_IO_IMPL_HPP

#include "./Display.hpp"
// TODO: Distribution plot
#include "./Print.hpp"
#include "./Read.hpp"
#include "./Spy.hpp"
#include "./Write.hpp"

namespace elem {

inline const char* QtImageFormat( FileFormat format )
{
    switch( format )
    {
    case BMP:  return "BMP";  break;
    case JPG:  return "JPG";  break;
    case JPEG: return "JPEG"; break;
    case PNG:  return "PNG";  break;
    case PPM:  return "PPM";  break;
    case XBM:  return "XBM";  break;
    case XPM:  return "XPM";  break;
    default: LogicError("Invalid image format"); return "N/A"; break;
    }
}

inline std::string FileExtension( FileFormat format )
{
    switch( format )
    {
    case ASCII:            return "txt";  break;
    case ASCII_MATLAB:     return "m";    break;
    case BINARY:           return "bin";  break;
    case BINARY_FLAT:      return "dat";  break;
    case BMP:              return "bmp";  break;
    case JPG:              return "jpg";  break;
    case JPEG:             return "jpeg"; break;
    case MATRIX_MARKET:    return "mm";   break;
    case PNG:              return "png";  break;
    case PPM:              return "ppm";  break;
    case XBM:              return "xbm";  break;
    case XPM:              return "xpm";  break;
    default: LogicError("Format not found"); return "N/A"; break;
    }
}

inline FileFormat FormatFromExtension( const std::string ext )
{
    bool foundFormat = false;
    FileFormat format = BINARY;
    for( int j=1; j<FileFormat_MAX; ++j )
    {
        format = static_cast<FileFormat>(j);
        if( FileExtension(format) == ext )
        {
            foundFormat = true;
            break;
        }
    }
    if( !foundFormat )
        RuntimeError("Did not detect file format");
    return format;
}

inline FileFormat DetectFormat( const std::string filename )
{
    const std::string ext = filename.substr(filename.find_last_of(".")+1);
    return FormatFromExtension( ext );
}

inline std::ifstream::pos_type
FileSize( std::ifstream& file )
{
    auto pos = file.tellg();
    file.seekg( 0, std::ifstream::end );
    auto numBytes = file.tellg();
    file.seekg( pos );
    return numBytes;
}

} // namespace elem

#endif // ifndef ELEM_IO_IMPL_HPP
