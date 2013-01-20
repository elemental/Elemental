/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_MEMORY_IMPL_HPP
#define CORE_MEMORY_IMPL_HPP

namespace elem {

template<typename G>
inline 
Memory<G>::Memory()
: size_(0), buffer_(NULL)
{ }

template<typename G>
inline 
Memory<G>::Memory( std::size_t size )
: size_(size), buffer_(new G[size])
{ }

template<typename G>
inline 
Memory<G>::~Memory()
{ delete[] buffer_; }

template<typename G>
inline G* 
Memory<G>::Buffer() const
{ return buffer_; }

template<typename G>
inline std::size_t 
Memory<G>::Size() const
{ return size_; }

template<typename G>
inline void 
Memory<G>::Require( std::size_t size )
{
    if( size > size_ )
    {
        delete[] buffer_;
#ifndef RELEASE
        try {
#endif
        buffer_ = new G[size];
#ifndef RELEASE
        } 
        catch( std::bad_alloc& exception )
        {
            std::cerr << "Failed to allocate " << size*sizeof(G) 
                      << " bytes" << std::endl;
            throw exception;
        }
#endif
        size_ = size;
    }
}

template<typename G>
inline void 
Memory<G>::Release()
{
#ifndef POOL_MEMORY
    this->Empty();
#endif
}

template<typename G>
inline void 
Memory<G>::Empty()
{
    delete[] buffer_;
    size_ = 0;
    buffer_ = 0;
}

} // namespace elem

#endif // ifndef CORE_MEMORY_IMPL_HPP
