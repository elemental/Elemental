/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_MEMORY_HPP
#define ELEMENTAL_MEMORY_HPP 1

namespace elemental {

template<typename T>
class Memory
{
    size_t _size;
    T*     _buffer;
public:
    Memory();
    Memory( size_t size );
    ~Memory();

    T*     Buffer() const;
    size_t Size()   const;

    void   Require( size_t size );
    void   Release();
};

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

template<typename T>
inline
elemental::Memory<T>::Memory()
: _size(0), _buffer(NULL)
{ }

template<typename T>
inline
elemental::Memory<T>::Memory( size_t size )
: _size(size), _buffer(new T[size])
{ }

template<typename T>
inline
elemental::Memory<T>::~Memory()
{ delete[] _buffer; }

template<typename T>
inline T*
elemental::Memory<T>::Buffer() const
{ return _buffer; }

template<typename T>
inline size_t
elemental::Memory<T>::Size() const
{ return _size; }

template<typename T>
inline void
elemental::Memory<T>::Require
( size_t size )
{
    if( size > _size )
    {
        delete[] _buffer;
        _buffer = new T[size];
        _size = size;
    }
}

template<typename T>
inline void
elemental::Memory<T>::Release()
{
#ifndef POOL_MEMORY
    delete[] _buffer;
    _size = 0;
    _buffer = 0;
#endif
}

#endif /* ELEMENTAL_MEMORY_HPP */

