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
#ifndef ELEMENTAL_LAPACK_INTERNAL_HPP
#define ELEMENTAL_LAPACK_INTERNAL_HPP 1

#include "elemental/lapack.hpp"

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

namespace elemental {
namespace lapack {
namespace internal {

//----------------------------------------------------------------------------//
// Local LAPACK                                                               //
//----------------------------------------------------------------------------//

template<typename F>
void
LocalChol
( Shape shape, DistMatrix<F,Star,Star>& A );

template<typename F>
void
LocalHegst
( Side side, Shape shape, 
  DistMatrix<F,Star,Star>& A, const DistMatrix<F,Star,Star>& B );

template<typename F>
void
LocalTrinv
( Shape shape, Diagonal diagonal, DistMatrix<F,Star,Star>& A );

//----------------------------------------------------------------------------//
// Chol helpers                                                               //
//----------------------------------------------------------------------------//

template<typename F>
void
CholLVar2( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholLVar2Naive( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholLVar3( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholLVar3Naive( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholLVar3Square( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholUVar2( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholUVar2Naive( DistMatrix<F,MC,MR>& A );
 
template<typename F>
void
CholUVar3( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholUVar3Naive( DistMatrix<F,MC,MR>& A );

template<typename F>
void
CholUVar3Square( DistMatrix<F,MC,MR>& A );
            
//----------------------------------------------------------------------------//
// GaussElim                                                                  //
//----------------------------------------------------------------------------//
            
template<typename F>
void
ReduceToRowEchelon
( DistMatrix<F,MC,MR>& A, DistMatrix<F,MC,MR>& B );

//----------------------------------------------------------------------------//
// Hegst                                                                      //
//----------------------------------------------------------------------------//

template<typename F>
void
HegstLLVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstLLVar2
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

// HegstLLVar3 would redundantly compute too much data

template<typename F>
void
HegstLLVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstLLVar5
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstLUVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstLUVar2
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

// HegstLUVar3 would redundantly compute too much data

template<typename F>
void
HegstLUVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstLUVar5
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstRLVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstRLVar2
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstRLVar3
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstRLVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstRLVar5
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L );

template<typename F>
void
HegstRUVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstRUVar2
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstRUVar3
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstRUVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

template<typename F>
void
HegstRUVar5
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U );

//----------------------------------------------------------------------------//
// LU                                                                         //
//----------------------------------------------------------------------------//

template<typename F>
void
ApplyRowPivots
(       DistMatrix<F,MC,MR>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage,
        int pivotOffset=0 );

void
ComposePivots
( const DistMatrix<int,Star,Star>& p,
        std::vector<int>& image,
        std::vector<int>& preimage,
        int pivotOffset = 0 );

template<typename F>
void
CreatePivotOp();

template<>
void
CreatePivotOp<float>();

template<>
void
CreatePivotOp<double>();

#ifndef WITHOUT_COMPLEX
template<>
void
CreatePivotOp<scomplex>();

template<>
void
CreatePivotOp<dcomplex>();
#endif

template<typename T>
void
DestroyPivotOp();

template<>
void
DestroyPivotOp<float>();

template<>
void
DestroyPivotOp<double>();

#ifndef WITHOUT_COMPLEX
template<>
void
DestroyPivotOp<scomplex>();

template<>
void
DestroyPivotOp<dcomplex>();
#endif

template<typename F>
void
LU( DistMatrix<F,MC,MR>& A, DistMatrix<int,Star,Star>& p );

template<typename F>
void
PanelLU
( DistMatrix<F,Star,Star>& A, 
  DistMatrix<F,MC,  Star>& B, 
  DistMatrix<int,Star,Star>& p, 
  int pivotOffset=0 );

template<typename F>
MPI_Op
PivotOp();

template<>
MPI_Op
PivotOp<float>();

template<>
MPI_Op
PivotOp<double>();

#ifndef WITHOUT_COMPLEX
template<>
MPI_Op
PivotOp<scomplex>();

template<>
MPI_Op
PivotOp<dcomplex>();
#endif
            
template<typename F>
void
PivotFunc
( void* inData, void* outData, 
  int* length, MPI_Datatype* datatype );

//----------------------------------------------------------------------------//
// QR                                                                         //
//----------------------------------------------------------------------------//

template<typename R>
void
PanelQR( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
PanelQR( DistMatrix<std::complex<R>,MC,MR  >& A,
         DistMatrix<std::complex<R>,MD,Star>& t );
#endif

//----------------------------------------------------------------------------//
// Reflector                                                                  //
//----------------------------------------------------------------------------//
template<typename R>
R 
ColReflector( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template<typename R>
std::complex<R>
ColReflector
( DistMatrix<std::complex<R>,MC,MR>& chi, 
  DistMatrix<std::complex<R>,MC,MR>& x );
#endif

template<typename R>
R
RowReflector( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template<typename R>
std::complex<R>
RowReflector
( DistMatrix<std::complex<R>,MC,MR>& chi,
  DistMatrix<std::complex<R>,MC,MR>& x );
#endif

//----------------------------------------------------------------------------//
// Tridiag                                                                    //
//----------------------------------------------------------------------------//

// Controls for which tridiagonalization method to use. 
//   TRIDIAG_NORMAL: 
//     Pipelined algorithm for nonsquare grids
//   TRIDIAG_SQUARE: 
//     Pipelined algorithm for square grids, drops down to a subset of processes
//     if necessary by redistributing the matrix.
//   TRIDIAG_DEFAULT:
//     Uses the TRIDIAG_NORMAL algorithm unless we are already on a square grid,
//     in which case the TRIDIAG_SQUARE algorithm is used.
//          
enum TridiagApproach { TRIDIAG_NORMAL, TRIDIAG_SQUARE, TRIDIAG_DEFAULT };
void SetTridiagApproach( TridiagApproach approach );

enum GridOrder { ROW_MAJOR, COL_MAJOR };
void SetTridiagSquareGridOrder( GridOrder order );

template<typename R>
void
PanelTridiagL
( DistMatrix<R,MC,MR  >& A, 
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,Star>& APan_MC_Star,
  DistMatrix<R,MR,Star>& APan_MR_Star,
  DistMatrix<R,MC,Star>& W_MC_Star,
  DistMatrix<R,MR,Star>& W_MR_Star );
template<typename R>
void
PanelTridiagU
( DistMatrix<R,MC,MR  >& A, 
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,Star>& APan_MC_Star,
  DistMatrix<R,MR,Star>& APan_MR_Star,
  DistMatrix<R,MC,Star>& W_MC_Star,
  DistMatrix<R,MR,Star>& W_MR_Star );
template<typename R>
void
PanelTridiagLSquare
( DistMatrix<R,MC,MR  >& A, 
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,Star>& APan_MC_Star,
  DistMatrix<R,MR,Star>& APan_MR_Star,
  DistMatrix<R,MC,Star>& W_MC_Star,
  DistMatrix<R,MR,Star>& W_MR_Star );
template<typename R>
void
PanelTridiagUSquare
( DistMatrix<R,MC,MR  >& A, 
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,Star>& APan_MC_Star,
  DistMatrix<R,MR,Star>& APan_MR_Star,
  DistMatrix<R,MC,Star>& W_MC_Star,
  DistMatrix<R,MR,Star>& W_MR_Star );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
PanelTridiagL
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,Star>& t,
  DistMatrix<std::complex<R>,MC,Star>& APan_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& APan_MR_Star,
  DistMatrix<std::complex<R>,MC,Star>& W_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& W_MR_Star );
template<typename R>
void
PanelTridiagU
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,Star>& t,
  DistMatrix<std::complex<R>,MC,Star>& APan_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& APan_MR_Star,
  DistMatrix<std::complex<R>,MC,Star>& W_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& W_MR_Star );

template<typename R>
void
PanelTridiagLSquare
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,Star>& t,
  DistMatrix<std::complex<R>,MC,Star>& APan_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& APan_MR_Star,
  DistMatrix<std::complex<R>,MC,Star>& W_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& W_MR_Star );
template<typename R>
void
PanelTridiagUSquare
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,Star>& t,
  DistMatrix<std::complex<R>,MC,Star>& APan_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& APan_MR_Star,
  DistMatrix<std::complex<R>,MC,Star>& W_MC_Star,
  DistMatrix<std::complex<R>,MR,Star>& W_MR_Star );
#endif
 
template<typename R>
void
TridiagL( DistMatrix<R,MC,MR>& A );
template<typename R>
void
TridiagU( DistMatrix<R,MC,MR>& A );

template<typename R>
void
TridiagLSquare( DistMatrix<R,MC,MR>& A );
template<typename R>
void
TridiagUSquare( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
TridiagL
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );
template<typename R>
void
TridiagU
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );

template<typename R>
void
TridiagLSquare
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );
template<typename R>
void
TridiagUSquare
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );
#endif

//----------------------------------------------------------------------------//
// Trinv                                                                      //
//----------------------------------------------------------------------------//

template<typename F>
void
TrinvVar3
( Shape shape, Diagonal diagonal, DistMatrix<F,MC,MR>& A  );

template<typename F>
void
TrinvLVar3
( Diagonal diagonal, DistMatrix<F,MC,MR>& L );

template<typename F>
void
TrinvUVar3
( Diagonal diagonal, DistMatrix<F,MC,MR>& U );

//----------------------------------------------------------------------------//
// UT Transform                                                               //
//----------------------------------------------------------------------------//

template<typename R>
void
UTLLN
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTLLN
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTLLC
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTLLC
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTLUN
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTLUN
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTLUC
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTLUC
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTRLN
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTRLN
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTRLC
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTRLC
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTRUN
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTRUN
( int offset,
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

template<typename R>
void
UTRUC
( int offset, 
  const DistMatrix<R,MC,MR>& H, 
        DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
UTRUC
( int offset, 
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,Star>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A );
#endif

//----------------------------------------------------------------------------//
// LAPACK Utility Functions                                                   //
//----------------------------------------------------------------------------//
template<typename F>
double
CholGFlops( int m, double seconds );

template<typename F>
double
HegstGFlops( int m, double seconds );

template<typename F>
double
LUGFlops( int m, double seconds );

template<typename F>
double
QRGFlops( int m, int n, double seconds );

template<typename F>
double
TridiagGFlops( int m, double seconds );

template<typename F>
double
TrinvGFlops( int m, double seconds );

template<typename F>
double
UTGFlops( int m, double seconds );

} // internal
} // lapack
} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

namespace elemental {
namespace lapack {
namespace internal {

//
// Local LAPACK
//

template<typename F>
inline void
LocalChol
( Shape shape, DistMatrix<F,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalChol");
#endif
    Chol( shape, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalHegst
( Side side, Shape shape,
  DistMatrix<F,Star,Star>& A, const DistMatrix<F,Star,Star>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalHegst");
#endif
    Hegst( side, shape, A.LocalMatrix(), B.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalTrinv
( Shape shape, Diagonal diagonal, DistMatrix<F,Star,Star>& A )
{ 
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalTrinv");
#endif
    Trinv( shape, diagonal, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// GFlop helpers
//

template<>
inline double
CholGFlops<float>
( int m, double seconds )
{ return (1./3.*m*m*m)/(1.e9*seconds); }
            
template<>
inline double
CholGFlops<double>
( int m, double seconds )
{ return CholGFlops<float>(m,seconds); }
            
#ifndef WITHOUT_COMPLEX
template<>
inline double
CholGFlops<scomplex>
( int m, double seconds )
{ return 4.*CholGFlops<float>(m,seconds); }
            
template<>
inline double
CholGFlops<dcomplex>
( int m, double seconds )
{ return 4.*CholGFlops<float>(m,seconds); }
#endif

template<>
inline double
HegstGFlops<float>
( int m, double seconds )
{ return (1.*m*m*m)/(1.e9*seconds); }

template<>
inline double
HegstGFlops<double>
( int m, double seconds )
{ return HegstGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
HegstGFlops<scomplex>
( int m, double seconds )
{ return 4.*HegstGFlops<float>(m,seconds); }

template<>
inline double
HegstGFlops<dcomplex>
( int m, double seconds )
{ return 4.*HegstGFlops<float>(m,seconds); }
#endif

template<>
inline double
LUGFlops<float>
( int m, double seconds )
{ return (2./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
LUGFlops<double>
( int m, double seconds )
{ return LUGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
LUGFlops<scomplex>
( int m, double seconds )
{ return 4.*LUGFlops<float>(m,seconds); }

template<>
inline double
LUGFlops<dcomplex>
( int m, double seconds )
{ return 4.*LUGFlops<float>(m,seconds); }
#endif

template<>
inline double
QRGFlops<float>
( int m, int n, double seconds )
{ return (2.*m*n*n-2./3.*n*n*n)/(1.e9*seconds); }

template<>
inline double
QRGFlops<double>
( int m, int n, double seconds )
{ return QRGFlops<float>(m,n,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
QRGFlops<scomplex>
( int m, int n, double seconds )
{ return 4.*QRGFlops<float>(m,n,seconds); }

template<>
inline double
QRGFlops<dcomplex>
( int m, int n, double seconds )
{ return 4.*QRGFlops<float>(m,n,seconds); }
#endif

template<>
inline double
TridiagGFlops<float>
( int m, double seconds )
{ return (4./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
TridiagGFlops<double>
( int m, double seconds )
{ return TridiagGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
TridiagGFlops<scomplex>
( int m, double seconds )
{ return 4.*TridiagGFlops<float>(m,seconds); }

template<>
inline double
TridiagGFlops<dcomplex>
( int m, double seconds )
{ return 4.*TridiagGFlops<float>(m,seconds); }
#endif

template<>
inline double
TrinvGFlops<float>
( int m, double seconds )
{ return (1./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
TrinvGFlops<double>
( int m, double seconds )
{ return TrinvGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
TrinvGFlops<scomplex>
( int m, double seconds )
{ return 4.*TrinvGFlops<float>(m,seconds); }

template<>
inline double
TrinvGFlops<dcomplex>
( int m, double seconds )
{ return 4.*TrinvGFlops<float>(m,seconds); }
#endif

template<>
inline double
UTGFlops<float>
( int m, double seconds )
{ return (2.*m*m*m)/(1.e9*seconds); }

template<>
inline double
UTGFlops<double>
( int m, double seconds )
{ return UTGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
UTGFlops<scomplex>
( int m, double seconds )
{ return 4.*UTGFlops<float>(m,seconds); }

template<>
inline double
UTGFlops<dcomplex>
( int m, double seconds )
{ return 4.*UTGFlops<float>(m,seconds); }
#endif

} // internal
} // lapack
} // elemental

#endif /* ELEMENTAL_LAPACK_INTERNAL_HPP */

