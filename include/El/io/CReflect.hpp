/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IO_CREFLECT_C_HPP
#define EL_IO_CREFLECT_C_HPP

namespace El {

// Input/Output
// ------------
inline ElColorMap CReflect( ColorMap map )
{ return static_cast<ElColorMap>(map); }
inline ColorMap CReflect( ElColorMap map )
{ return static_cast<ColorMap>(map); }

} // namespace El

#endif // ifndef EL_IO_CREFLECT_C_HPP
