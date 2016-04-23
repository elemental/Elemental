/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_PROPS_CREFLECT_C_HPP
#define EL_LAPACK_PROPS_CREFLECT_C_HPP

namespace El {

inline ElNormType CReflect( NormType type )
{ return static_cast<ElNormType>(type); }
inline NormType CReflect( ElNormType type )
{ return static_cast<NormType>(type); }

} // namespace El

#endif // ifndef EL_LAPACK_PROPS_CREFLECT_C_HPP
