/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_BLAS_H
#define ELEMENTAL_BLAS_H 1

#include <cmath>
#include "ElementalDistMatrix.h"
#include "ElementalPartitioning.h"
#include "wrappers/BLAS.h"

namespace Elemental 
{
    namespace BLAS 
    {
        //--------------------------------------------------------------------//
        // Local BLAS: Level 0 (extensions)                                   //
        //--------------------------------------------------------------------//

        template<typename R>
        R
        Abs( const R& alpha );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R 
        Abs( const std::complex<R>& alpha );
#endif

        template<typename R>
        R
        FastAbs( const R& alpha );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        FastAbs( const std::complex<R>& alpha );
#endif
        
        template<typename R>
        R 
        Conj( const R& alpha );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        std::complex<R>
        Conj( const std::complex<R>& alpha );
#endif

        template<typename R>
        R
        Imag( const R& alpha );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        Imag( const std::complex<R>& alpha );
#endif

        template<typename R>
        R
        Real( const R& alpha );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        Real( const std::complex<R>& alpha );
#endif

        //--------------------------------------------------------------------//
        // Local BLAS: Level 1                                                //
        //--------------------------------------------------------------------//

        // AXPY: Y := Alpha X Plus Y 
        template<typename T>
        void
        Axpy( const T alpha, const Matrix<T>& X, Matrix<T>& Y );

        // COPY: Copy
        template<typename T>
        void
        Copy( const Matrix<T>& X, Matrix<T>& Y );

        // DOT: alpha := conj(x)^T * y
        // 
        // Though the standard BLAS interface only defines DOT for real 
        // datatypes, it is naturally generalized to an inner product over the
        // complex field. Recall that the conjugate symmetry of inner products 
        // requires that (x,y) = conj(y,x), so that (x,x) = conj( (x,x) ) => 
        // (x,x) is real. This requires that we choose (x,x) = conj(x)^T * x.
        template<typename T>
        T
        Dot( const Matrix<T>& x, const Matrix<T>& y );

        // DOTC: alpha := conj(x)^T * y
        //
        // This is the sister routine to DOT; while DOT is originally defined 
        // only over the reals, DOTC was defined only over the complex field. 
        // They are each others' extensions, and so, to us, they are 
        // identical.
        template<typename T>
        T
        Dotc( const Matrix<T>& x, const Matrix<T>& y );

        // DOTU: alpha := x^T * y
        //
        // Standard BLAS defines DOTU for complex datatypes, but the operation
        // is perfectly valid over the reals (clearly), so we extend it.
        template<typename T>
        T
        Dotu( const Matrix<T>& x, const Matrix<T>& y );

        // NRM2: NoRM 2 (Euclidean norm)
        template<typename R>
        R
        Nrm2( const Matrix<R>& x ); 

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        Nrm2( const Matrix< std::complex<R> >& x );
#endif

        // SCAL: SCALe X by alpha
        template<typename T>
        void
        Scal( const T alpha, Matrix<T>& X );
        
        //--------------------------------------------------------------------//
        // Local BLAS: Level 1 (extensions)                                   //
        //--------------------------------------------------------------------//

        // CONJ: CONJugate in-place
        //
        // There are two implementations because, for real datatypes, Conj is
        // a no-op. Partial specialization of function templates is not allowed,
        // so we must have two declarations.
        template<typename R>
        void
        Conj( Matrix<R>& A );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        void
        Conj( Matrix< std::complex<R> >& A );
#endif

        // CONJ: CONJugated copy
        template<typename T>
        void
        Conj( const Matrix<T>& A, Matrix<T>& B );

        // CONJTRANS: CONJugated Transposed copy
        template<typename T>
        void
        ConjTrans( const Matrix<T>& A, Matrix<T>& B );

        // TRANS: TRANSposed copy
        template<typename T>
        void
        Trans( const Matrix<T>& A, Matrix<T>& B );

        //--------------------------------------------------------------------//
        // Local BLAS: Level 2                                                //
        //--------------------------------------------------------------------//

        // GEMV: GEneral Matrix-Vector multiply
        template<typename T>
        void
        Gemv( const Orientation orientation,
              const T alpha, const Matrix<T>& A, const Matrix<T>& x,
              const T beta,        Matrix<T>& y                     );

        // GER: GEneral Rank-one update
        //
        // For complex datatypes it routes to Gerc, as x (tensor product) y 
        // is x * conj(y)^T. That is, the dual of y is its conjugate transpose
        // thanks to the Riesz map. Thus our generalized Ger is equivalent to
        // our generalized Gerc.
        template<typename T>
        void
        Ger( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                  Matrix<T>& A                     );

        // GERC: GEneral Rank-one Conjugated update
        template<typename T>
        void
        Gerc( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        // GERU: GEneral Rank-one Unconjugated update
        template<typename T>
        void
        Geru( const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        // HEMV: HErmitian Matrix-Vector multiply
        template<typename T>
        void
        Hemv( const Shape shape,
              const T alpha, const Matrix<T>& A, const Matrix<T>& x,
              const T beta,        Matrix<T>& y                     );

        // HER: HErmitian Rank-one update
        template<typename T>
        void
        Her( const Shape shape,
             const T alpha, const Matrix<T>& x, Matrix<T>& A );

        // HER2: HErmitian Rank-2 update
        template<typename T>
        void
        Her2( const Shape shape,
              const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        // SYMV: SYmmetric Matrix-Vector multiply
        template<typename T>
        void
        Symv( const Shape shape,
              const T alpha, const Matrix<T>& A, const Matrix<T>& x,
              const T beta,        Matrix<T>& y                     );

        // SYR: SYmmetric Rank-one update
        template<typename T>
        void
        Syr( const Shape shape,
             const T alpha, const Matrix<T>& x, Matrix<T>& A );

        // SYR2: SYmmetric Rank-2 update
        template<typename T>
        void
        Syr2( const Shape shape,
              const T alpha, const Matrix<T>& x, const Matrix<T>& y,
                                   Matrix<T>& A                     );

        // TRMV: TRiangular Matrix-Vector multiply
        template<typename T>
        void
        Trmv
        ( const Shape shape, 
          const Orientation orientation, 
          const Diagonal diagonal, 
          const Matrix<T>& A, 
                Matrix<T>& x            );

        // TRSV: TRiangular Solve with a Vector
        template<typename T>
        void
        Trsv
        ( const Shape shape, 
          const Orientation orientation,
          const Diagonal diagonal,
          const Matrix<T>& A, 
                Matrix<T>& x            );

        //--------------------------------------------------------------------//
        // Local BLAS: Level 3                                                //
        //--------------------------------------------------------------------//

        // GEMM: GEneral Matrix-Matrix multiplication
        template<typename T>
        void
        Gemm
        ( const Orientation orientationOfA, const Orientation orientationOfB,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta,        Matrix<T>& C                                  );

        // HEMM: HErmitian Matrix-Matrix multiply
        template<typename T>
        void
        Hemm
        ( const Side side, const Shape shape,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta,        Matrix<T>& C                     );

        // HER2K: HErmitian Rank-2K update
        template<typename T>
        void
        Her2k
        ( const Shape shape, const Orientation orientation,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta, Matrix<T>& C );

        // HERK: HErmitian Rank-K update
        template<typename T>
        void
        Herk
        ( const Shape shape, const Orientation orientation,
          const T alpha, const Matrix<T>& A, const T beta, Matrix<T>& C );

        // SYMM: SYmmetric Matrix-Matrix multiply
        template<typename T>
        void
        Symm
        ( const Side side, const Shape shape,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta,        Matrix<T>& C                     );

        // SYR2K: SYmmetric Rank-2K update
        template<typename T>
        void
        Syr2k
        ( const Shape shape, const Orientation orientation,
          const T alpha, const Matrix<T>& A, const Matrix<T>& B,
          const T beta, Matrix<T>& C );

        // SYRK: SYmmetric Rank-K update
        template<typename T>
        void
        Syrk
        ( const Shape shape, const Orientation orientation,
          const T alpha, const Matrix<T>& A, const T beta, Matrix<T>& C );

        // TRMM: TRiangular Matrix-Matrix multiplication
        template<typename T>
        void
        Trmm
        ( const Side side,
          const Shape shape,
          const Orientation orientation,
          const Diagonal diagonal,
          const T alpha, const Matrix<T>& A, Matrix<T>& B );

        // TRSM: TRiangular Solve with Multiple right-hand sides
        template<typename T>
        void
        Trsm
        ( const Side side, 
          const Shape shape,
          const Orientation orientation,
          const Diagonal diagonal,
          const T alpha, const Matrix<T>& A, Matrix<T>& B ); 
        
        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 1                                          //
        //--------------------------------------------------------------------//

        // AXPY: Y := Alpha X Plus Y 
        template<typename T, Distribution U, Distribution V>
        void
        Axpy( const T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y );

        // COPY: Copy
        //
        // In our case, it is just a wrapper around the '=' operator for those
        // that prefer BLAS/PLAPACK syntax.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        // DOT: alpha := conj(x)^T * y
        // 
        // Though the standard BLAS interface only defines DOT for real 
        // datatypes, it is naturally generalized to an inner product over the
        // complex field. Recall that the conjugate symmetry of inner products 
        // requires that (x,y) = conj(y,x), so that (x,x) = conj( (x,x) ) => 
        // (x,x) is real. This requires that we choose (x,x) = conj(x)^T * x.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dot( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // DOTC: alpha := conj(x)^T * y
        //
        // This is the sister routine to DOT; while DOT is originally defined 
        // only over the reals, DOTC was defined only over the complex field. 
        // They are each others' extensions, and so, to us, they are 
        // identical.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dotc( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // DOTU: alpha := x^T * y
        //
        // Standard BLAS defines DOTU for complex datatypes, but the operation
        // is perfectly valid over the reals (clearly), so we extend it.
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        T
        Dotu( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y );

        // NRM2: NoRM 2 (Euclidean norm)
        template<typename R>
        R
        Nrm2( const DistMatrix<R,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
        template<typename R>
        R
        Nrm2( const DistMatrix< std::complex<R>, MC, MR >& x );
#endif

        // SCAL: SCALe by a constant
        template<typename T, Distribution U, Distribution V>
        void
        Scal
        ( const T alpha, DistMatrix<T,U,V>& A );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 1 (extensions)                             //
        //--------------------------------------------------------------------//

        // CONJ: CONJugate in-place
        template<typename T, Distribution U, Distribution V>
        void
        Conj( DistMatrix<T,U,V>& A );

        // CONJ: CONJugated copy
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        Conj( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        // CONJTRANS: CONJugated Transposed copy
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        ConjTrans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        // TRANS: TRANSposed copy
        template<typename T, Distribution U, Distribution V,
                             Distribution W, Distribution Z >
        void
        Trans( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 2                                          //
        //--------------------------------------------------------------------//

        // GEMV: GEneral Matrix-Vector multiplication
        template<typename T>
        void
        Gemv
        ( const Orientation orientationOfA,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& x,
          const T beta,        DistMatrix<T,MC,MR>& y );

        // GER: GEneral Rank-one update
        //
        // For complex datatypes it routes to Gerc, as x (tensor product) y 
        // is x * conj(y)^T. That is, the dual of y is its conjugate transpose
        // thanks to the Riesz map. Thus our generalized Ger is equivalent to
        // our generalized Gerc.
        template<typename T>
        void
        Ger
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );
        
        // GERC: GEneral Rank-one Conjugated update
        //
        // Since the extension of Ger to complex datatypes 
        template<typename T>
        void
        Gerc
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );
        
        // GERU: GEneral Rank-one Unconjugated update
        template<typename T>
        void
        Geru
        ( const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );

        // HEMV: HErmitian Matrix-Vector multiplication
        template<typename T>
        void
        Hemv
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& x,
          const T beta,        DistMatrix<T,MC,MR>& y );

        // HER: HErmitian Rank-one update
        template<typename T>
        void
        Her
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& x,
                               DistMatrix<T,MC,MR>& A );
        
        // HER: HErmitian Rank-2 update
        template<typename T>
        void
        Her2
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );

        // SYMV: SYmmetric Matrix-Vector multiplication
        template<typename T>
        void
        Symv
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& x,
          const T beta,        DistMatrix<T,MC,MR>& y );

        // SYR: SYmmetric Rank-one update
        template<typename T>
        void
        Syr
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& x,
                               DistMatrix<T,MC,MR>& A );

        // SYR2: SYmmetric Rank-2 update
        template<typename T>
        void
        Syr2
        ( const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& x,
                         const DistMatrix<T,MC,MR>& y,
                               DistMatrix<T,MC,MR>& A );

        // TRMV: TRiangular Matrix-Vector multiply
        template<typename T>
        void
        Trmv
        ( const Shape shape, 
          const Orientation orientation, 
          const Diagonal diagonal,
          const DistMatrix<T,MC,MR>& A, 
                DistMatrix<T,MC,MR>& x  );
        
        // TRSV: TRiangular Solve with a Vector
        template<typename T>
        void
        Trsv
        ( const Shape shape, 
          const Orientation orientation, 
          const Diagonal diagonal,
          const DistMatrix<T,MC,MR>& A, 
                DistMatrix<T,MC,MR>& x  );

        //--------------------------------------------------------------------//
        // Distributed BLAS: Level 3                                          //
        //--------------------------------------------------------------------//

        // GEMM: GEneral Matrix-Matrix multiplication
        template<typename T>
        void
        Gemm
        ( const Orientation orientationOfA, 
          const Orientation orientationOfB,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C );

        // HEMM: HErmitian Matrix-Matrix multiply
        template<typename T>
        void
        Hemm
        ( const Side side, const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C );

        // HER2K: HErmitian Rank-2K Update
        template<typename T>
        void
        Her2k
        ( const Shape shape, const Orientation orientation,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C     );

        // HERK: HErmitian Rank-K Update
        template<typename T>
        void
        Herk
        ( const Shape shape, const Orientation orientation,
          const T alpha, const DistMatrix<T,MC,MR>& A,
          const T beta,        DistMatrix<T,MC,MR>& C      );

        // SYMM: SYmmetric Matrix-Matrix multiply
        template<typename T>
        void
        Symm
        ( const Side side, const Shape shape,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C );

        // SYR2K: SYmmetric Rank-2K Update
        template<typename T>
        void
        Syr2k
        ( const Shape shape, const Orientation orientation,
          const T alpha, const DistMatrix<T,MC,MR>& A,
                         const DistMatrix<T,MC,MR>& B,
          const T beta,        DistMatrix<T,MC,MR>& C      );

        // SYRK: SYmmetric Rank-K Update
        template<typename T>
        void
        Syrk
        ( const Shape shape, const Orientation orientation,
          const T alpha, const DistMatrix<T,MC,MR>& A,
          const T beta,        DistMatrix<T,MC,MR>& C      );

        // TRMM: TRiangular Matrix-Matrix multiplication
        template<typename T>
        void
        Trmm
        ( const Side side, 
          const Shape shape,
          const Orientation orientation, 
          const Diagonal diagonal,
          const T alpha, 
          const DistMatrix<T,MC,MR>& A,
                DistMatrix<T,MC,MR>& B  );

        // TRSM: TRiangular Solve with Multiplie right-hand sides
        template<typename T>
        void
        Trsm
        ( const Side side,
          const Shape shape,
          const Orientation orientation,
          const Diagonal diagonal,
          const T alpha, 
          const DistMatrix<T,MC,MR>& A,
                DistMatrix<T,MC,MR>& B  );
    }
}

/*----------------------------------------------------------------------------*/

//----------------------------------------------------------------------------//
// Local BLAS: Level 0 (extensions)                                           //
//----------------------------------------------------------------------------//

template<typename R>
inline R
Elemental::BLAS::Abs
( const R& alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::BLAS::Abs
( const std::complex<R>& alpha )
{ return std::abs( alpha ); }
#endif

template<typename R>
inline R
Elemental::BLAS::FastAbs
( const R& alpha )
{ return fabs(alpha); }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::BLAS::FastAbs
( const std::complex<R>& alpha )
{ return fabs( std::real(alpha) ) + fabs( std::imag(alpha) ); }
#endif

template<typename R>
inline R
Elemental::BLAS::Conj
( const R& alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline std::complex<R>
Elemental::BLAS::Conj
( const std::complex<R>& alpha )
{ return std::conj( alpha ); }
#endif

template<typename R>
inline R
Elemental::BLAS::Imag
( const R& alpha )
{ return 0; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::BLAS::Imag
( const std::complex<R>& alpha )
{ return std::imag( alpha ); }
#endif

template<typename R>
inline R
Elemental::BLAS::Real
( const R& alpha )
{ return alpha; }

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::BLAS::Real
( const std::complex<R>& alpha )
{ return std::real( alpha ); }
#endif

//----------------------------------------------------------------------------//
// Local BLAS: Level 1                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Axpy
( const T alpha, const Matrix<T>& X, Matrix<T>& Y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Axpy");
    if( X.Height() != Y.Height() || X.Width() != Y.Width() )
    {
        std::cerr << "Nonconformal Axpy." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    if( X.Width() <= X.Height() )
    {
        for( int j=0; j<X.Width(); ++j )
        {
            Elemental::wrappers::BLAS::Axpy
            ( X.Height(), alpha, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
        }
    }
    else
    {
        for( int i=0; i<X.Height(); ++i )
        {
            Elemental::wrappers::BLAS::Axpy
            ( X.Width(), alpha, X.LockedBuffer(i,0), X.LDim(),
                                Y.Buffer(i,0),       Y.LDim() );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Copy
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline T
Elemental::BLAS::Dot
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Dot");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        std::cerr << "Expected vector inputs." << std::endl;        
        DumpCallStack();
        throw std::exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        std::cerr << "x and y must be the same length." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dot
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::BLAS::Dot
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Elemental::BLAS::Dotc
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Dotc");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        std::cerr << "Expected vector inputs." << std::endl;        
        DumpCallStack();
        throw std::exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        std::cerr << "x and y must be the same length." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotc
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::BLAS::Dotc
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename T>
inline T
Elemental::BLAS::Dotu
( const Matrix<T>& x, const Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Dotu");
    if( (x.Height() != 1 && x.Width() != 1) ||
        (y.Height() != 1 && y.Width() != 1)   )
    {
        std::cerr << "Expected vector inputs." << std::endl;        
        DumpCallStack();
        throw std::exception();
    }
    int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    int yLength = ( y.Width() == 1 ? y.Height() : y.Width() );
    if( xLength != yLength )
    {
        std::cerr << "x and y must be the same length." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    T dotProduct;
    if( x.Width() == 1 && y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), 1 );
    }
    else if( x.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotu
                     ( x.Height(), x.LockedBuffer(), 1,
                                   y.LockedBuffer(), y.LDim() );
    }
    else if( y.Width() == 1 )
    {
        dotProduct = wrappers::BLAS::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), 1        );
    }
    else
    {
        dotProduct = wrappers::BLAS::Dotu
                     ( x.Width(), x.LockedBuffer(), x.LDim(),
                                  y.LockedBuffer(), y.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return dotProduct;
}

template<typename R>
inline R
Elemental::BLAS::Nrm2
( const Matrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("BLAS::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
    {
        std::cerr << "Expected vector input." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = Elemental::wrappers::BLAS::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = Elemental::wrappers::BLAS::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

#ifndef WITHOUT_COMPLEX
template<typename R>
inline R
Elemental::BLAS::Nrm2
( const Matrix< std::complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("BLAS::Nrm2");
    if( x.Height() != 1 && x.Width() != 1 )
    {
        std::cerr << "Expected vector input." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    R norm;
    if( x.Width() == 1 )
    {
        norm = Elemental::wrappers::BLAS::Nrm2
               ( x.Height(), x.LockedBuffer(), 1 );
    }
    else
    {
        norm = Elemental::wrappers::BLAS::Nrm2
               ( x.Width(), x.LockedBuffer(), x.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}
#endif

template<typename T>
inline void
Elemental::BLAS::Scal
( const T alpha, Matrix<T>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Scal");
#endif
    if( alpha != (T)1 )
    {
        if( X.Width() <= X.Height() )
        {
            for( int j=0; j<X.Width(); ++j )
            {
                Elemental::wrappers::BLAS::Scal
                ( X.Height(), alpha, X.Buffer(0,j), 1 );
            }
        }
        else
        {
            for( int i=0; i<X.Height(); ++i )
            {
                Elemental::wrappers::BLAS::Scal
                ( X.Width(), alpha, X.Buffer(i,0), X.LDim() );
            }
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 1 (extensions)                                           //
//----------------------------------------------------------------------------//

// Default case is for real datatypes
template<typename R>
inline void
Elemental::BLAS::Conj
( Matrix<R>& A )
{ }

#ifndef WITHOUT_COMPLEX
// Specialization is to complex datatypes
template<typename R>
inline void
Elemental::BLAS::Conj
( Matrix< std::complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj (in-place)");
#endif
    const int m = A.Height();
    const int n = A.Width();
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            A(i,j) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template<typename T>
inline void
Elemental::BLAS::Conj
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(i,j) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::ConjTrans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::ConjTrans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = BLAS::Conj( A(i,j) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trans
( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trans");
#endif
    const int m = A.Height();
    const int n = A.Width();
    B.ResizeTo( n, m );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            B(j,i) = A(i,j);
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 2                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Gemv
( const Orientation orientation,
  const T alpha, const Matrix<T>& A, const Matrix<T>& x,
  const T beta,        Matrix<T>& y                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemv");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( orientation == Normal )
    {
        if( A.Height() != yLength || A.Width() != xLength )
        {
            std::cerr << "A must conform with x and y:" << 
            std::endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            std::endl << "  x ~ " << x.Height() << " x " << x.Width() <<
            std::endl << "  y ~ " << y.Height() << " x " << y.Width() << 
            std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Width() != yLength || A.Height() != xLength )
        {
            std::cerr << "A must conform with x and y:" << 
            std::endl << "  A ~ " << A.Height() << " x " << A.Width() <<
            std::endl << "  x ~ " << x.Height() << " x " << x.Width() <<
            std::endl << "  y ~ " << y.Height() << " x " << y.Width() << 
            std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char transChar = OrientationToChar( orientation );
    const int m = A.Height();
    const int n = A.Width();
    const int k = ( transChar == 'N' ? n : m );
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    if( k != 0 )
    {
        Elemental::wrappers::BLAS::Gemv
        ( transChar, m, n, 
          alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
          beta,  y.Buffer(), incy                                   );
    }
    else
    {
        Elemental::BLAS::Scal( beta, y );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Ger
( const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Ger");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::cerr << "Nonconformal Ger: " << std::endl;
        std::cerr << "  x ~ " << x.Height() << " x " << x.Width() << std::endl;
        std::cerr << "  y ~ " << y.Height() << " x " << y.Width() << std::endl;
        std::cerr << "  A ~ " << A.Height() << " x " << A.Width() << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Ger
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Gerc
( const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gerc");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::cerr << "Nonconformal Gerc." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Gerc
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Geru
( const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Geru");
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Width() )
    {
        std::cerr << "Nonconformal Geru." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const int m = A.Height(); 
    const int n = A.Width();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Geru
    ( m, n, alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy, 
                   A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Hemv
( const Shape shape,
  const T alpha, const Matrix<T>& A, const Matrix<T>& x,
  const T beta,        Matrix<T>& y                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Hemv");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        std::cerr << "A must conform with x and y." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Hemv
    ( uploChar, m,
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx,
      beta,  y.Buffer(), incy                                   );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Her
( const Shape shape,
  const T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Her");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( x.Width() != 1 && x.Height() != 1 )
    {
        std::cerr << "x must be a vector." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
    {
        std::cerr << "x must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    Elemental::wrappers::BLAS::Her
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Her2
( const Shape shape,
  const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Her2");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1)   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
    {
        std::cerr << "x and y must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Her2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Symv
( const Shape shape,
  const T alpha, const Matrix<T>& A, const Matrix<T>& x,
  const T beta,        Matrix<T>& y                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Symv");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( ( x.Height() != 1 && x.Width() != 1 ) ||
        ( y.Height() != 1 && y.Width() != 1 )   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( A.Height() != xLength || A.Height() != yLength )
    {
        std::cerr << "A must conform with x and y." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Symv
    ( uploChar, m, 
      alpha, A.LockedBuffer(), A.LDim(), x.LockedBuffer(), incx, 
      beta,  y.Buffer(), incy                                   );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Syr
( const Shape shape,
  const T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syr");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( x.Width() != 1 && x.Height() != 1 )
    {
        std::cerr << "x must be a vector." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
    {
        std::cerr << "x must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    Elemental::wrappers::BLAS::Syr
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, A.Buffer(), A.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Syr2
( const Shape shape,
  const T alpha, const Matrix<T>& x, const Matrix<T>& y, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syr2");
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( (x.Width() != 1 && x.Height() != 1) || 
        (y.Width() != 1 && y.Height() != 1)   )
    {
        std::cerr << "x and y must be vectors." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    const int yLength = ( y.Width()==1 ? y.Height() : y.Width() );
    if( xLength != A.Height() || yLength != A.Height() )
    {
        std::cerr << "x and y must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const int incy = ( y.Width()==1 ? 1 : y.LDim() );
    Elemental::wrappers::BLAS::Syr2
    ( uploChar, m,
      alpha, x.LockedBuffer(), incx, y.LockedBuffer(), incy,
             A.Buffer(), A.LDim()                           );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trmv
( const Shape shape, const Orientation orientation, const Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x                                          )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trmv");
    if( x.Height() != 1 && x.Width() != 1 )
    {
        std::cerr << "x must be a vector." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack(); 
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
    {
        std::cerr << "x must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    Elemental::wrappers::BLAS::Trmv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trsv
( const Shape shape, const Orientation orientation, const Diagonal diagonal,
  const Matrix<T>& A, Matrix<T>& x                                          )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trsv");
    if( x.Height() != 1 && x.Width() != 1 )
    {
        std::cerr << "x must be a vector." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( A.Height() != A.Width() )
    {
        std::cerr << "A must be square." << std::endl;
        DumpCallStack(); 
        throw std::exception();
    }
    const int xLength = ( x.Width()==1 ? x.Height() : x.Width() );
    if( xLength != A.Height() )
    {
        std::cerr << "x must conform with A." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    const int m = A.Height();
    const int incx = ( x.Width()==1 ? 1 : x.LDim() );
    Elemental::wrappers::BLAS::Trsv
    ( uploChar, transChar, diagChar, m,
      A.LockedBuffer(), A.LDim(), x.Buffer(), incx );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Local BLAS: Level 3                                                        //
//----------------------------------------------------------------------------//

template<typename T>
inline void
Elemental::BLAS::Gemm
( const Orientation orientationOfA, const Orientation orientationOfB,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B,
  const T beta,        Matrix<T>& C                                  )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gemm");
    if( orientationOfA == Normal && orientationOfB == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Width()  != B.Height()    )
        {
            std::cerr << "Nonconformal GemmNN." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientationOfA == Normal )
    {
        if( A.Height() != C.Height() ||
            B.Height() != C.Width()  ||
            A.Width()  != B.Width()     )
        {
            std::cerr << "Nonconformal GemmN(T/C)." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientationOfB == Normal )
    {
        if( A.Width()  != C.Height() ||
            B.Width()  != C.Width()  ||
            A.Height() != B.Height()    )
        {
            std::cerr << "Nonconformal Gemm(T/C)N." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Width()  != C.Height() ||
            B.Height() != C.Width()  ||
            A.Height() != B.Width()     )
        {
            std::cerr << "Nonconformal Gemm(T/C)(T/C)." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char transA = OrientationToChar( orientationOfA );
    const char transB = OrientationToChar( orientationOfB );
    const int m = C.Height();
    const int n = C.Width();
    const int k = ( orientationOfA == Normal ? A.Width() : A.Height() );
    if( k != 0 )
    {
        Elemental::wrappers::BLAS::Gemm
        ( transA, transB, m, n, k, 
          alpha, A.LockedBuffer(), A.LDim(), B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        Elemental::BLAS::Scal( beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Hemm
( const Side side, const Shape shape,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B, 
  const T beta,        Matrix<T>& C                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Hemm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    wrappers::BLAS::Hemm( sideChar, shapeChar, C.Height(), C.Width(),
                          alpha, A.LockedBuffer(), A.LDim(), 
                                 B.LockedBuffer(), B.LDim(),
                          beta,  C.Buffer(),       C.LDim()          );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Her2k
( const Shape shape, const Orientation orientation,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B,
  const T beta,        Matrix<T>& C                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Her2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
        {
            std::cerr << "Nonconformal Her2k." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
        {
            std::cerr << "Nonconformal Her2k." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        std::cerr << "Her2k only accepts Normal and ConjugateTranpose options." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::BLAS::Her2k( uplo, trans, C.Height(), k, 
                           alpha, A.LockedBuffer(), A.LDim(), 
                                  B.LockedBuffer(), B.LDim(),
                           beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Herk
( const Shape shape, const Orientation orientation,
  const T alpha, const Matrix<T>& A, const T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Herk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
        {
            std::cerr << "Nonconformal Herk." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientation == ConjugateTranspose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
        {
            std::cerr << "Nonconformal Herk." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        std::cerr << "Herk only accepts Normal and ConjugateTranpose options." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::BLAS::Herk( uplo, trans, C.Height(), k, 
                          alpha, A.LockedBuffer(), A.LDim(), 
                          beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Symm
( const Side side, const Shape shape,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B, 
  const T beta,        Matrix<T>& C                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Symm");
#endif
    const char sideChar = SideToChar( side );
    const char shapeChar = ShapeToChar( shape );
    wrappers::BLAS::Symm( sideChar, shapeChar, C.Height(), C.Width(),
                          alpha, A.LockedBuffer(), A.LDim(), 
                                 B.LockedBuffer(), B.LDim(),
                          beta,  C.Buffer(),       C.LDim()          );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Syr2k
( const Shape shape, const Orientation orientation,
  const T alpha, const Matrix<T>& A, const Matrix<T>& B,
  const T beta,        Matrix<T>& C                     )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syr2k");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
        {
            std::cerr << "Nonconformal Syr2k." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
        {
            std::cerr << "Nonconformal Syr2k." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        std::cerr << "Syr2k only accepts Normal and Tranpose options." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::BLAS::Syr2k( uplo, trans, C.Height(), k, 
                           alpha, A.LockedBuffer(), A.LDim(), 
                                  B.LockedBuffer(), B.LDim(),
                           beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Syrk
( const Shape shape, const Orientation orientation,
  const T alpha, const Matrix<T>& A, const T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("BLAS::Syrk");
    if( orientation == Normal )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
        {
            std::cerr << "Nonconformal Syrk." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else if( orientation == Transpose )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
        {
            std::cerr << "Nonconformal Syrk." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        std::cerr << "Syrk only accepts Normal and Tranpose options." 
                  << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    const char uplo = ShapeToChar( shape );
    const char trans = OrientationToChar( orientation );
    const int k = ( orientation == Normal ? A.Width() : A.Height() );
    wrappers::BLAS::Syrk( uplo, trans, C.Height(), k, 
                          alpha, A.LockedBuffer(), A.LDim(), 
                          beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trmm
( const Side side, 
  const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trmm");
    if( A.Height() != A.Width() )
    {
        std::cerr << "Triangular matrix must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( side == Left )
    {
        if( A.Height() != B.Height() )
        {
            std::cerr << "Nonconformal Trmm." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Height() != B.Width() )
        {
            std::cerr << "Nonconformal Trmm." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    Elemental::wrappers::BLAS::Trmm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Elemental::BLAS::Trsm
( const Side side, 
  const Shape shape,
  const Orientation orientation,
  const Diagonal diagonal,
  const T alpha, const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trsm");
    if( A.Height() != A.Width() )
    {
        std::cerr << "Triangular matrix must be square." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
    if( side == Left )
    {
        if( A.Height() != B.Height() )
        {
            std::cerr << "Nonconformal Trsm." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
    else
    {
        if( A.Height() != B.Width() )
        {
            std::cerr << "Nonconformal Trsm." << std::endl;
            DumpCallStack();
            throw std::exception();
        }
    }
#endif
    const char sideChar = SideToChar( side );
    const char uploChar = ShapeToChar( shape );
    const char transChar = OrientationToChar( orientation );
    const char diagChar = DiagonalToChar( diagonal );
    Elemental::wrappers::BLAS::Trsm
    ( sideChar, uploChar, transChar, diagChar, B.Height(), B.Width(),
      alpha, A.LockedBuffer(), A.LDim(), B.Buffer(), B.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 1                                                  //
//----------------------------------------------------------------------------//

template<typename T, Elemental::Distribution U, Elemental::Distribution V>
inline void
Elemental::BLAS::Axpy
( const T alpha, const DistMatrix<T,U,V>& X, DistMatrix<T,U,V>& Y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Axpy");
    const Grid& grid = X.GetGrid();
    if( X.GetGrid() != Y.GetGrid() )
    {
        if( grid.VCRank() == 0 ) 
        {
            std::cerr << "X and Y must be distributed over the same grid."
                      << std::endl;
        }
        DumpCallStack();
        throw std::exception();
    }
    if( X.ColAlignment() != Y.ColAlignment() ||
        X.RowAlignment() != Y.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
            std::cerr << "Axpy requires X and Y be aligned." << std::endl;
        DumpCallStack();
        throw std::exception();
    }
#endif
    BLAS::Axpy( alpha, X.LockedLocalMatrix(), Y.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::Copy
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Copy");
#endif
    B = A;
#ifndef RELEASE
    PopCallStack();
#endif
}

// Our extended Dotc is equivalent to our extended Dot, 
// but we are burdened with consistency
template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline T
Elemental::BLAS::Dotc
( const DistMatrix<T,U,V>& x, const DistMatrix<T,W,Z>& y )
{
#ifndef RELEASE
    PushCallStack("BLAS::Dotc");
#endif
    T globalDot = BLAS::Dot( x, y );
#ifndef RELEASE
    PopCallStack();
#endif
    return globalDot;
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V>
inline void
Elemental::BLAS::Scal
( const T alpha, DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Scal");
#endif
    BLAS::Scal( alpha, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}


//----------------------------------------------------------------------------//
// Distributed BLAS: Level 1 (extensions)                                     //
//----------------------------------------------------------------------------//

template<typename T, Elemental::Distribution U, Elemental::Distribution V>
inline void
Elemental::BLAS::Conj
( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj (in-place)");
#endif
    BLAS::Conj( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::Conj
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Conj");
#endif
    B = A;
    BLAS::Conj( B ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::ConjTrans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::ConjTrans");
#endif
    DistMatrix<T,Z,W> C( B.GetGrid() );
    if( B.ConstrainedColDist() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowDist() )
        C.AlignColsWith( B );

    C = A;

    if( !B.ConstrainedColDist() )
        B.AlignColsWith( C );
    if( !B.ConstrainedRowDist() )
        B.AlignRowsWith( C );

    B.ResizeTo( A.Width(), A.Height() );
    BLAS::ConjTrans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T, Elemental::Distribution U, Elemental::Distribution V,
                     Elemental::Distribution W, Elemental::Distribution Z >
inline void
Elemental::BLAS::Trans
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    PushCallStack("BLAS::Trans");
#endif
    DistMatrix<T,Z,W> C( B.GetGrid() );
    if( B.ConstrainedColDist() )
        C.AlignRowsWith( B );
    if( B.ConstrainedRowDist() )
        C.AlignColsWith( B );

    C = A;

    if( !B.ConstrainedColDist() )
        B.AlignColsWith( C );
    if( !B.ConstrainedRowDist() )
        B.AlignRowsWith( C );

    B.ResizeTo( A.Width(), A.Height() );
    BLAS::Trans( C.LockedLocalMatrix(), B.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 2                                                  //
//----------------------------------------------------------------------------//

// Our extended Ger and Gerc are equivalent
template<typename T>
inline void
Elemental::BLAS::Gerc
( const T alpha, const DistMatrix<T,MC,MR>& x,
                 const DistMatrix<T,MC,MR>& y,
                       DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("BLAS::Gerc");
#endif
    BLAS::Ger( alpha, x, y, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

#endif /* ELEMENTAL_BLAS_H */

