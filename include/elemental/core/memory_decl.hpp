/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

template<typename G>
class Memory
{
    std::size_t size_;
    G*     buffer_;
public:
    Memory();
    Memory( std::size_t size );
    ~Memory();

    G*          Buffer() const;
    std::size_t Size()   const;

    void   Require( std::size_t size );
    void   Release();
    void   Empty();
};

} // namespace elem

