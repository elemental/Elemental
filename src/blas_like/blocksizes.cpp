/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like.hpp>
#include <stack>

namespace {
using namespace El;

std::stack<Int> blocksizeStack;

template<typename T>
struct LocalSymvBlocksizeHelper { static Int value; };
template<typename T>
Int LocalSymvBlocksizeHelper<T>::value = 64;

template<typename T>
struct LocalTrrkBlocksizeHelper { static Int value; };
template<typename T>
Int LocalTrrkBlocksizeHelper<T>::value = 64;

template<typename T>
struct LocalTrr2kBlocksizeHelper { static Int value; };
template<typename T>
Int LocalTrr2kBlocksizeHelper<T>::value = 64;

}

namespace El {

Int Blocksize()
{
    EL_DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to extract blocksize from empty stack");
    )
    return ::blocksizeStack.top();
}

void SetBlocksize( Int blocksize )
{
    EL_DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to set blocksize at top of empty stack");
    )
    ::blocksizeStack.top() = blocksize;
}

void PushBlocksizeStack( Int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void PopBlocksizeStack()
{
    EL_DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to pop an empty blocksize stack");
    )
    ::blocksizeStack.pop();
}

void EmptyBlocksizeStack()
{
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
}

template<typename T>
void SetLocalSymvBlocksize( Int blocksize )
{ LocalSymvBlocksizeHelper<T>::value = blocksize; }

template<typename T>
Int LocalSymvBlocksize()
{ return LocalSymvBlocksizeHelper<T>::value; }

template<typename T>
void SetLocalTrrkBlocksize( Int blocksize )
{ LocalTrrkBlocksizeHelper<T>::value = blocksize; }

template<typename T>
Int LocalTrrkBlocksize()
{ return LocalTrrkBlocksizeHelper<T>::value; }

template<typename T>
void SetLocalTrr2kBlocksize( Int blocksize )
{ LocalTrr2kBlocksizeHelper<T>::value = blocksize; }

template<typename T>
Int LocalTrr2kBlocksize()
{ return LocalTrr2kBlocksizeHelper<T>::value; }

#define PROTO(T) \
  template void SetLocalSymvBlocksize<T>( Int blocksize ); \
  template Int LocalSymvBlocksize<T>(); \
  template void SetLocalTrrkBlocksize<T>( Int blocksize ); \
  template Int LocalTrrkBlocksize<T>(); \
  template void SetLocalTrr2kBlocksize<T>( Int blocksize ); \
  template Int LocalTrr2kBlocksize<T>();

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
