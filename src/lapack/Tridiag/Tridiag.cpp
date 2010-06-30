/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

template<typename R>
void
elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    if( shape == Lower )
        lapack::internal::TridiagL( A, d, e );
    else
        lapack::internal::TridiagU( A, d, e );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<R,         MD,Star>& d,
  DistMatrix<R,         MD,Star>& e )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    if( shape == Lower )
        lapack::internal::TridiagL( A, d, e );
    else
        lapack::internal::TridiagU( A, d, e );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MD,Star>& d,
  DistMatrix<float,MD,Star>& e );

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MD,Star>& d,
  DistMatrix<double,MD,Star>& e );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<float,   MD,Star>& d,
  DistMatrix<float,   MD,Star>& e );

template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<double,  MD,Star>& d,
  DistMatrix<double,  MD,Star>& e );
#endif

