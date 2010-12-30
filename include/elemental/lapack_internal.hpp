/*
   Copyright (c) 2009-2010, Jack Poulson
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

namespace elemental {
namespace lapack {
namespace internal {

//----------------------------------------------------------------------------//
// Local LAPACK                                                               //
//----------------------------------------------------------------------------//

template<typename T>
void
LocalChol
( Shape shape, DistMatrix<T,Star,Star>& A );

template<typename T>
void
LocalHegst
( Side side, Shape shape, 
  DistMatrix<T,Star,Star>& A, const DistMatrix<T,Star,Star>& B );

template<typename T>
void
LocalTrinv
( Shape shape, Diagonal diagonal, DistMatrix<T,Star,Star>& A );

//----------------------------------------------------------------------------//
// Chol helpers                                                               //
//----------------------------------------------------------------------------//

template<typename T>
void
CholVar2
( Shape shape, DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholVar2Naive
( Shape shape, DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholVar3
( Shape shape, DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholVar3Naive
( Shape shape, DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholL( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholU( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholLVar2( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholLVar2Naive( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholLVar3( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholLVar3Naive( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholUVar2( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholUVar2Naive( DistMatrix<T,MC,MR>& A );
 
template<typename T>
void
CholUVar3( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholUVar3Naive( DistMatrix<T,MC,MR>& A );
            
//----------------------------------------------------------------------------//
// GaussElim                                                                  //
//----------------------------------------------------------------------------//
            
template<typename T>
void
ReduceToRowEchelon
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

//----------------------------------------------------------------------------//
// Hegst                                                                      //
//----------------------------------------------------------------------------//

template<typename T>
void
HegstOld
( Side side, Shape shape,
  DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B );

template<typename T>
void
HegstNaive
( Side side, Shape shape, 
  DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B );

template<typename T>
void
HegstLL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstLLNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstLLOld
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstLU
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstLUNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstLUOld
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstRL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstRLNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstRLOld
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstRU
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstRUNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstRUOld
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

//----------------------------------------------------------------------------//
// LU                                                                         //
//----------------------------------------------------------------------------//

template<typename T>
void
ApplyRowPivots
(       DistMatrix<T,MC,MR>& A, 
  const std::vector<int>& image,
  const std::vector<int>& preimage,
        int pivotOffset=0 );

void
ComposePivots
( const DistMatrix<int,Star,Star>& p,
        std::vector<int>& image,
        std::vector<int>& preimage,
        int pivotOffset = 0 );

template<typename T>
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

template<typename T>
void
LU( DistMatrix<T,MC,MR>& A, DistMatrix<int,Star,Star>& p );

template<typename T>
void
PanelLU
( DistMatrix<T,Star,Star>& A, 
  DistMatrix<T,VC,  Star>& B, 
  DistMatrix<int,Star,Star>& p, 
  int pivotOffset=0 );

template<typename T>
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
            
template<typename T>
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
PanelTridiagLSquare
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
PanelTridiagLSquare
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
PanelTridiagU
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e );
#ifndef WITHOUT_COMPLEX
template<typename R>
void
PanelTridiagU
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
#endif

template<typename R>
void
TridiagL( DistMatrix<R,MC,MR>& A );
template<typename R>
void
TridiagLSquare( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
TridiagL
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );
template<typename R>
void
TridiagLSquare
( DistMatrix<std::complex<R>,MC,  MR  >& A, 
  DistMatrix<std::complex<R>,Star,Star>& t );
#endif

template<typename R>
void
TridiagU( DistMatrix<R,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
TridiagU
( DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,Star,Star>& t );
#endif

//----------------------------------------------------------------------------//
// Trinv                                                                      //
//----------------------------------------------------------------------------//

template<typename T>
void
TrinvVar3
( Shape shape, Diagonal diagonal, DistMatrix<T,MC,MR>& A  );

template<typename T>
void
TrinvL
( Diagonal diagonal, DistMatrix<T,MC,MR>& L );

template<typename T>
void
TrinvU
( Diagonal diagonal, DistMatrix<T,MC,MR>& U );

template<typename T>
void
TrinvLVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& L );

template<typename T>
void
TrinvUVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& U );

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
template<typename T>
double
CholGFlops( int m, double seconds );

template<typename T>
double
HegstGFlops( int m, double seconds );

template<typename T>
double
LUGFlops( int m, double seconds );

template<typename T>
double
QRGFlops( int m, int n, double seconds );

template<typename T>
double
TridiagGFlops( int m, double seconds );

template<typename T>
double
TrinvGFlops( int m, double seconds );

template<typename T>
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

template<typename T>
inline void
LocalChol
( Shape shape, DistMatrix<T,Star,Star>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalChol");
#endif
    Chol( shape, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalHegst
( Side side, Shape shape,
  DistMatrix<T,Star,Star>& A, const DistMatrix<T,Star,Star>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalHegst");
#endif
    Hegst( side, shape, A.LocalMatrix(), B.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalTrinv
( Shape shape, Diagonal diagonal, DistMatrix<T,Star,Star>& A )
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

