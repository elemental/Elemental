/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
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
( bool bothOnLeft, Shape shape, 
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
CholVar3
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
CholLVar3( DistMatrix<T,MC,MR>& A );

template<typename T>
void
CholUVar2( DistMatrix<T,MC,MR>& A );
 
template<typename T>
void
CholUVar3( DistMatrix<T,MC,MR>& A );
            
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
HegstFalseL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstFalseU
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U );

template<typename T>
void
HegstTrueL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L );

template<typename T>
void
HegstTrueU
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
// Tridiag                                                                    //
//----------------------------------------------------------------------------//
template<typename R>
R 
Reflector( DistMatrix<R,MC,MR>& chi, DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template<typename R>
std::complex<R>
Reflector
( DistMatrix<std::complex<R>,MC,MR>& chi, 
  DistMatrix<std::complex<R>,MC,MR>& x );
#endif

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

template<typename R>
void
PanelTridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
PanelTridiagL
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
#endif
 
template<typename R>
void
PanelTridiagU
( DistMatrix<R,MC,  MR  >& A,
  DistMatrix<R,MC,  MR  >& W,
  DistMatrix<R,MD,  Star>& e,
  DistMatrix<R,MD,  Star>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
PanelTridiagU
( DistMatrix<std::complex<R>,MC,  MR  >& A,
  DistMatrix<std::complex<R>,MC,  MR  >& W,
  DistMatrix<R,              MD,  Star>& e,
  DistMatrix<std::complex<R>,MD,  Star>& t );
#endif

template<typename R>
void
TridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
TridiagL
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<R,              MD,Star>& d,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
#endif

template<typename R>
void
TridiagU
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template<typename R>
void
TridiagU
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<R,              MD,Star>& d,
  DistMatrix<R,              MD,Star>& e,
  DistMatrix<std::complex<R>,MD,Star>& t );
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
TridiagGFlops( int m, double seconds );

template<typename T>
double
TrinvGFlops( int m, double seconds );

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
( bool bothOnLeft, Shape shape,
  DistMatrix<T,Star,Star>& A, const DistMatrix<T,Star,Star>& B )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::LocalHegst");
#endif
    Hegst( bothOnLeft, shape, A.LocalMatrix(), B.LockedLocalMatrix() );
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

} // internal
} // lapack
} // elemental

#endif /* ELEMENTAL_LAPACK_INTERNAL_HPP */

