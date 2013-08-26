/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_MEMORY_DECL_HPP
#define ELEM_CORE_MEMORY_DECL_HPP

namespace elem {

template<typename G>
class Memory
{
    std::size_t size_;
    G* buffer_;
public:
    Memory();
    Memory( std::size_t size );
    ~Memory();

#ifndef SWIG
    Memory( Memory<G>&& mem );
    Memory<G>& operator=( Memory<G>&& mem );
#endif
    void Swap( Memory<G>& mem );

    G* Buffer() const;
    std::size_t Size()   const;

    G* Require( std::size_t size );
    void Release();
    void Empty();
};

} // namespace elem

#endif // ifndef ELEM_CORE_MEMORY_DECL_HPP
