/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// Generate a sample from a uniform PDF over the (closed) unit ball about the 
// origin of the ring implied by the type T using the most natural metric.
template<typename T> T SampleUnitBall();

} // namespace elem
