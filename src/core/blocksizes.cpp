/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <stack>

// TODO: Consider moving this into another folder (e.g., in src/blas_like/)
//       and/or using a better datastructure for managing the blocksizes

namespace {
using namespace El;

std::stack<Int> blocksizeStack;

// Tuning parameters for basic routines
Int localSymvIntBlocksize = 64;
Int localSymvFloatBlocksize = 64;
Int localSymvDoubleBlocksize = 64;
Int localSymvComplexFloatBlocksize = 64;
Int localSymvComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localSymvDoubleDoubleBlocksize = 64;
Int localSymvQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localSymvQuadBlocksize = 64;
Int localSymvComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localSymvBigIntBlocksize = 64;
Int localSymvBigFloatBlocksize = 64;
#endif

Int localTrr2kIntBlocksize = 64;
Int localTrr2kFloatBlocksize = 64;
Int localTrr2kDoubleBlocksize = 64;
Int localTrr2kComplexFloatBlocksize = 64;
Int localTrr2kComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localTrr2kDoubleDoubleBlocksize = 64;
Int localTrr2kQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localTrr2kQuadBlocksize = 64;
Int localTrr2kComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localTrr2kBigIntBlocksize = 64;
Int localTrr2kBigFloatBlocksize = 64;
#endif

Int localTrrkIntBlocksize = 64;
Int localTrrkFloatBlocksize = 64;
Int localTrrkDoubleBlocksize = 64;
Int localTrrkComplexFloatBlocksize = 64;
Int localTrrkComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localTrrkDoubleDoubleBlocksize = 64;
Int localTrrkQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localTrrkQuadBlocksize = 64;
Int localTrrkComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localTrrkBigIntBlocksize = 64;
Int localTrrkBigFloatBlocksize = 64;
#endif

}

namespace El {

Int Blocksize()
{ 
    DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to extract blocksize from empty stack");
    )
    return ::blocksizeStack.top(); 
}

void SetBlocksize( Int blocksize )
{ 
    DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to set blocksize at top of empty stack");
    )
    ::blocksizeStack.top() = blocksize; 
}

void PushBlocksizeStack( Int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void PopBlocksizeStack()
{
    DEBUG_ONLY(
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

template<>
void SetLocalSymvBlocksize<Int>( Int blocksize )
{ ::localSymvIntBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<float>( Int blocksize )
{ ::localSymvFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<double>( Int blocksize )
{ ::localSymvDoubleBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<float>>( Int blocksize )
{ ::localSymvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<double>>( Int blocksize )
{ ::localSymvComplexDoubleBlocksize = blocksize; }

#ifdef EL_HAVE_QD
template<>
void SetLocalSymvBlocksize<DoubleDouble>( Int blocksize )
{ ::localSymvDoubleDoubleBlocksize = blocksize; }
template<>
void SetLocalSymvBlocksize<QuadDouble>( Int blocksize )
{ ::localSymvQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalSymvBlocksize<Quad>( Int blocksize )
{ ::localSymvQuadBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<Quad>>( Int blocksize )
{ ::localSymvComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalSymvBlocksize<BigInt>( Int blocksize )
{ ::localSymvBigIntBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<BigFloat>( Int blocksize )
{ ::localSymvBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalSymvBlocksize<Int>()
{ return ::localSymvIntBlocksize; }

template<>
Int LocalSymvBlocksize<float>()
{ return ::localSymvFloatBlocksize; }

template<>
Int LocalSymvBlocksize<double>()
{ return ::localSymvDoubleBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<float>>()
{ return ::localSymvComplexFloatBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<double>>()
{ return ::localSymvComplexDoubleBlocksize; }

#ifdef EL_HAVE_QD
template<>
Int LocalSymvBlocksize<DoubleDouble>()
{ return ::localSymvDoubleDoubleBlocksize; }
template<>
Int LocalSymvBlocksize<QuadDouble>()
{ return ::localSymvQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalSymvBlocksize<Quad>()
{ return ::localSymvQuadBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<Quad>>()
{ return ::localSymvComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalSymvBlocksize<BigInt>()
{ return ::localSymvBigIntBlocksize; }
template<>
Int LocalSymvBlocksize<BigFloat>()
{ return ::localSymvBigFloatBlocksize; }
#endif

template<>
void SetLocalTrr2kBlocksize<Int>( Int blocksize )
{ ::localTrr2kIntBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<float>( Int blocksize )
{ ::localTrr2kFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<double>( Int blocksize )
{ ::localTrr2kDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<float>>( Int blocksize )
{ ::localTrr2kComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<double>>( Int blocksize )
{ ::localTrr2kComplexDoubleBlocksize = blocksize; }

#ifdef EL_HAVE_QD
template<>
void SetLocalTrr2kBlocksize<DoubleDouble>( Int blocksize )
{ ::localTrr2kDoubleDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<QuadDouble>( Int blocksize )
{ ::localTrr2kQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalTrr2kBlocksize<Quad>( Int blocksize )
{ ::localTrr2kQuadBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<Quad>>( Int blocksize )
{ ::localTrr2kComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalTrr2kBlocksize<BigInt>( Int blocksize )
{ ::localTrr2kBigIntBlocksize = blocksize; }
template<>
void SetLocalTrr2kBlocksize<BigFloat>( Int blocksize )
{ ::localTrr2kBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalTrr2kBlocksize<Int>()
{ return ::localTrr2kIntBlocksize; }

template<>
Int LocalTrr2kBlocksize<float>()
{ return ::localTrr2kFloatBlocksize; }

template<>
Int LocalTrr2kBlocksize<double>()
{ return ::localTrr2kDoubleBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<float>>()
{ return ::localTrr2kComplexFloatBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<double>>()
{ return ::localTrr2kComplexDoubleBlocksize; }

#ifdef EL_HAVE_QD
template<>
Int LocalTrr2kBlocksize<DoubleDouble>()
{ return ::localTrr2kDoubleDoubleBlocksize; }

template<>
Int LocalTrr2kBlocksize<QuadDouble>()
{ return ::localTrr2kQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalTrr2kBlocksize<Quad>()
{ return ::localTrr2kQuadBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<Quad>>()
{ return ::localTrr2kComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalTrr2kBlocksize<BigInt>()
{ return ::localTrr2kBigIntBlocksize; }
template<>
Int LocalTrr2kBlocksize<BigFloat>()
{ return ::localTrr2kBigFloatBlocksize; }
#endif

template<>
void SetLocalTrrkBlocksize<Int>( Int blocksize )
{ ::localTrrkIntBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<float>( Int blocksize )
{ ::localTrrkFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<double>( Int blocksize )
{ ::localTrrkDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<float>>( Int blocksize )
{ ::localTrrkComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<double>>( Int blocksize )
{ ::localTrrkComplexDoubleBlocksize = blocksize; }

#ifdef EL_HAVE_QD
template<>
void SetLocalTrrkBlocksize<DoubleDouble>( Int blocksize )
{ ::localTrrkDoubleDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<QuadDouble>( Int blocksize )
{ ::localTrrkQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalTrrkBlocksize<Quad>( Int blocksize )
{ ::localTrrkQuadBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<Quad>>( Int blocksize )
{ ::localTrrkComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalTrrkBlocksize<BigInt>( Int blocksize )
{ ::localTrrkBigIntBlocksize = blocksize; }
template<>
void SetLocalTrrkBlocksize<BigFloat>( Int blocksize )
{ ::localTrrkBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalTrrkBlocksize<Int>()
{ return ::localTrrkIntBlocksize; }

template<>
Int LocalTrrkBlocksize<float>()
{ return ::localTrrkFloatBlocksize; }

template<>
Int LocalTrrkBlocksize<double>()
{ return ::localTrrkDoubleBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<float>>()
{ return ::localTrrkComplexFloatBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<double>>()
{ return ::localTrrkComplexDoubleBlocksize; }

#ifdef EL_HAVE_QD
template<>
Int LocalTrrkBlocksize<DoubleDouble>()
{ return ::localTrrkDoubleDoubleBlocksize; }

template<>
Int LocalTrrkBlocksize<QuadDouble>()
{ return ::localTrrkQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalTrrkBlocksize<Quad>()
{ return ::localTrrkQuadBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<Quad>>()
{ return ::localTrrkComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalTrrkBlocksize<BigInt>()
{ return ::localTrrkBigIntBlocksize; }
template<>
Int LocalTrrkBlocksize<BigFloat>()
{ return ::localTrrkBigFloatBlocksize; }
#endif

} // namespace El
