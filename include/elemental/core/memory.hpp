/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_MEMORY_HPP
#define ELEMENTAL_MEMORY_HPP 1

namespace elemental {

template<typename G>
class Memory
{
    size_t _size;
    G*     _buffer;
public:
    Memory();
    Memory( size_t size );
    ~Memory();

    G*     Buffer() const;
    size_t Size()   const;

    void   Require( size_t size );
    void   Release();
    void   Empty();
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename G>
inline Memory<G>::Memory()
: _size(0), _buffer(NULL)
{ }

template<typename G>
inline Memory<G>::Memory( size_t size )
: _size(size), _buffer(new G[size])
{ }

template<typename G>
inline Memory<G>::~Memory()
{ delete[] _buffer; }

template<typename G>
inline G* Memory<G>::Buffer() const
{ return _buffer; }

template<typename G>
inline size_t Memory<G>::Size() const
{ return _size; }

template<typename G>
inline void Memory<G>::Require( size_t size )
{
    if( size > _size )
    {
        delete[] _buffer;
        _buffer = new G[size];
        _size = size;
    }
}

template<typename G>
inline void Memory<G>::Release()
{
#ifndef POOL_MEMORY
    this->Empty();
#endif
}

template<typename G>
inline void Memory<G>::Empty()
{
    delete[] _buffer;
    _size = 0;
    _buffer = 0;
}

} // namespace elemental

#endif /* ELEMENTAL_MEMORY_HPP */

