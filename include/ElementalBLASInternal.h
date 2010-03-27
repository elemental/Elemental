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
#ifndef ELEMENTAL_BLAS_INTERNAL_H
#define ELEMENTAL_BLAS_INTERNAL_H 1

#include "ElementalBLAS.h"

namespace Elemental
{
    namespace BLAS
    {
        namespace Internal
        {
            //----------------------------------------------------------------//
            // Distributed BLAS: Level 1                                      //
            //----------------------------------------------------------------//
            
            // Pseudo-partial-specializations of BLAS::Dot
            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MC,MR>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MC,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,MR>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MR,MC>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MR,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,MC>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,VC,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,VC>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,VR,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,VR>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dot
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,Star>& y );

            // Pseudo-partial-specializations of BLAS::Dotu
            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MC,MR>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MC,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,MR>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MR,MC>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,MR,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,MC>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,VC,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,VC>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,VR,Star>& y );

            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,VR>& y );
            
            template<typename T, Distribution U, Distribution V>
            T
            Dotu
            ( const DistMatrix<T,U,V>& x, 
              const DistMatrix<T,Star,Star>& y );

            //----------------------------------------------------------------//
            // Distributed BLAS: Level 2                                      //
            //----------------------------------------------------------------//
            
            // Gemv where A is not transposed
            template<typename T>
            void
            GemvN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& x,
              const T beta,        DistMatrix<T,MC,MR>& y );

            // Gemv where A is (conjugate) transposed
            template<typename T>
            void
            GemvT
            ( const Orientation orientation,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& x,
              const T beta,        DistMatrix<T,MC,MR>& y );

            // This is for the case where x is a column vector.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^H x[MC,* ]
            template<typename T>
            void
            HemvColAccumulate
            ( const Shape shape,
              const T alpha, 
              const DistMatrix<T,MC,MR  >& A,
              const DistMatrix<T,MC,Star>& x_MC_Star,
              const DistMatrix<T,MR,Star>& x_MR_Star,
                    DistMatrix<T,MC,Star>& z_MC_Star,
                    DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a column vector and A is lower.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^H x[MC,* ]
            template<typename T>
            void
            HemvColAccumulateL
            ( const T alpha, const DistMatrix<T,MC,MR  >& A,
                             const DistMatrix<T,MC,Star>& x_MC_Star,
                             const DistMatrix<T,MR,Star>& x_MR_Star,
                                   DistMatrix<T,MC,Star>& z_MC_Star,
                                   DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a column vector and A is upper.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha triu(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trius(A)[MC,MR])^H x[MC,* ]
            template<typename T>
            void
            HemvColAccumulateU
            ( const T alpha, const DistMatrix<T,MC,MR  >& A,
                             const DistMatrix<T,MC,Star>& x_MC_Star,
                             const DistMatrix<T,MR,Star>& x_MR_Star,
                                   DistMatrix<T,MC,Star>& z_MC_Star,
                                   DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a row vector.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^H (x[* ,MC])^T
            template<typename T>
            void
            HemvRowAccumulate
            ( const Shape shape,
              const T alpha, 
              const DistMatrix<T,MC,  MR  >& A,
              const DistMatrix<T,Star,MC  >& x_Star_MC,
              const DistMatrix<T,Star,MR  >& x_Star_MR,
                    DistMatrix<T,MC,  Star>& z_MC_Star,
                    DistMatrix<T,MR,  Star>& z_MR_Star
            );

            // This is for the case where x is a row vector and A is lower.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^H (x[* ,MC])^T
            template<typename T>
            void
            HemvRowAccumulateL
            ( const T alpha, 
              const DistMatrix<T,MC,  MR  >& A,
              const DistMatrix<T,Star,MC  >& x_Star_MC,
              const DistMatrix<T,Star,MR  >& x_Star_MR,
                    DistMatrix<T,MC,  Star>& z_MC_Star,
                    DistMatrix<T,MR,  Star>& z_MR_Star
            );

            // This is for the case where x is a row vector and A is upper.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha triu(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trius(A)[MC,MR])^H (x[* ,MC])^T
            template<typename T>
            void
            HemvRowAccumulateU
            ( const T alpha, 
              const DistMatrix<T,MC,  MR  >& A,
              const DistMatrix<T,Star,MC  >& x_Star_MC,
              const DistMatrix<T,Star,MR  >& x_Star_MR,
                    DistMatrix<T,MC,  Star>& z_MC_Star,
                    DistMatrix<T,MR,  Star>& z_MR_Star
            );


            // This is for the case where x is a column vector.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^T x[MC,* ]
            template<typename T>
            void
            SymvColAccumulate
            ( const Shape shape,
              const T alpha, 
              const DistMatrix<T,MC,MR  >& A,
              const DistMatrix<T,MC,Star>& x_MC_Star,
              const DistMatrix<T,MR,Star>& x_MR_Star,
                    DistMatrix<T,MC,Star>& z_MC_Star,
                    DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a column vector and A is lower.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^T x[MC,* ]
            template<typename T>
            void
            SymvColAccumulateL
            ( const T alpha, 
              const DistMatrix<T,MC,MR  >& A,
              const DistMatrix<T,MC,Star>& x_MC_Star,
              const DistMatrix<T,MR,Star>& x_MR_Star,
                    DistMatrix<T,MC,Star>& z_MC_Star,
                    DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a column vector and A is upper.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha triu(A)[MC,MR] x[MR,* ]
            //     z[MR,* ] := alpha (trius(A)[MC,MR])^T x[MC,* ]
            template<typename T>
            void
            SymvColAccumulateU
            ( const T alpha, 
              const DistMatrix<T,MC,MR  >& A,
              const DistMatrix<T,MC,Star>& x_MC_Star,
              const DistMatrix<T,MR,Star>& x_MR_Star,
                    DistMatrix<T,MC,Star>& z_MC_Star,
                    DistMatrix<T,MR,Star>& z_MR_Star
            );

            // This is for the case where x is a row vector.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^T (x[* ,MC])^T
            template<typename T>
            void
            SymvRowAccumulate
            ( const Shape shape,
              const T alpha, 
              const DistMatrix<T,MC,  MR>& A,
              const DistMatrix<T,Star,MC>& x_Star_MC,
              const DistMatrix<T,Star,MR>& x_Star_MR,
                    DistMatrix<T,Star,MC>& z_Star_MC,
                    DistMatrix<T,Star,MR>& z_Star_MR
            );

            // This is for the case where x is a row vector and A is lower.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trils(A)[MC,MR])^T (x[* ,MC])^T
            template<typename T>
            void
            SymvRowAccumulateL
            ( const T alpha, 
              const DistMatrix<T,MC,  MR>& A,
              const DistMatrix<T,Star,MC>& x_Star_MC,
              const DistMatrix<T,Star,MR>& x_Star_MR,
                    DistMatrix<T,Star,MC>& z_Star_MC,
                    DistMatrix<T,Star,MR>& z_Star_MR
            );

            // This is for the case where x is a row vector and A is upper.
            //
            // Returns the unreduced components z[MC,* ] and z[MR,* ]:
            //     z[MC,* ] := alpha triu(A)[MC,MR] (x[* ,MR])^T
            //     z[MR,* ] := alpha (trius(A)[MC,MR])^T (x[* ,MC])^T
            template<typename T>
            void
            SymvRowAccumulateU
            ( const T alpha, 
              const DistMatrix<T,MC,  MR>& A,
              const DistMatrix<T,Star,MC>& x_Star_MC,
              const DistMatrix<T,Star,MR>& x_Star_MR,
                    DistMatrix<T,Star,MC>& z_Star_MC,
                    DistMatrix<T,Star,MR>& z_Star_MR
            );

            template<typename T>
            void
            TrsvLN
            ( const Diagonal diagonal,
              const DistMatrix<T,MC,MR>& L,
                    DistMatrix<T,MC,MR>& x );

            template<typename T>
            void
            TrsvLT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const DistMatrix<T,MC,MR>& L,
                    DistMatrix<T,MC,MR>& x );

            template<typename T>
            void
            TrsvUN
            ( const Diagonal diagonal,
              const DistMatrix<T,MC,MR>& U,
                    DistMatrix<T,MC,MR>& x );

            template<typename T>
            void
            TrsvUT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const DistMatrix<T,MC,MR>& U,
                    DistMatrix<T,MC,MR>& x  );

            //----------------------------------------------------------------//
            // Distributed BLAS: Level 3                                      //
            //----------------------------------------------------------------//

            // Gemm where we avoid redistributing A.
            template<typename T>
            void
            GemmA
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta,
                    DistMatrix<T,MC,MR>& C    );

            // Gemm where we avoid redistributing B.
            template<typename T>
            void
            GemmB
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta, 
                    DistMatrix<T,MC,MR>& C     );

            // Gemm where we avoid redistributing C.
            template<typename T>
            void
            GemmC
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, 
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta,
                    DistMatrix<T,MC,MR>& C     );

            // Gemm for panel-panel dot products.
            template<typename T>
            void
            GemmDot
            ( const Orientation orientationOfA,
              const Orientation orientationOfB,
              const T alpha,
              const DistMatrix<T,MC,MR>& A,
              const DistMatrix<T,MC,MR>& B,
              const T beta,
                    DistMatrix<T,MC,MR>& C     );

            // Normal Normal Gemm.
            template<typename T>
            void
            GemmNN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing A.
            template<typename T>
            void
            GemmNNA
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing B.
            template<typename T>
            void
            GemmNNB
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm where we avoid redistributing C.
            template<typename T>
            void
            GemmNNC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal Normal Gemm for panel-panel dot product
            template<typename T>
            void
            GemmNNDot
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal (Conjugate)Transpose Gemm.
            template<typename T>
            void
            GemmNT
            ( const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal (Conjugate)Transpose Gemm where we avoid redistributing A.
            template<typename T>
            void
            GemmNTA
            ( const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal (Conjugate)Transpose Gemm where we avoid redistributing B.
            template<typename T>
            void
            GemmNTB
            ( const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Normal (Conjugate)Transpose Gemm where we avoid redistributing C.
            template<typename T>
            void
            GemmNTC
            ( const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );
 
            // Normal (Conjugate)Transpose Gemm for panel-panel dot product
            template<typename T>
            void
            GemmNTDot
            ( const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );
 
            // (Conjugate)Transpose Normal Gemm.
            template<typename T>
            void
            GemmTN
            ( const Orientation orientationOfA,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose Normal Gemm where we avoid redistributing A.
            template<typename T>
            void
            GemmTNA
            ( const Orientation orientationOfA,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose Normal Gemm where we avoid redistributing B.
            template<typename T>
            void
            GemmTNB
            ( const Orientation orientationOfA,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose Normal Gemm where we avoid redistributing C.
            template<typename T>
            void
            GemmTNC
            ( const Orientation orientationOfA,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose Normal Gemm for panel-panel dot product
            template<typename T>
            void
            GemmTNDot
            ( const Orientation orientationOfA,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose (Conjugate)Transpose Gemm.
            template<typename T>
            void
            GemmTT
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
            // redistributing A.
            template<typename T>
            void
            GemmTTA
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
            // redistributing B.
            template<typename T>
            void
            GemmTTB
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
            // redistributing C.
            template<typename T>
            void
            GemmTTC
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // (Conjugate)Transpose (Conjugate)Transpose Gemm for panel-panel
            // dot product
            template<typename T>
            void
            GemmTTDot
            ( const Orientation orientationOfA, 
              const Orientation orientationOfB,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Hemm
            template<typename T>
            void
            Hemm
            ( const Side side, const Shape shape,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Lower Hemm
            template<typename T>
            void
            HemmLL
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Lower Hemm where we avoid redistributing C
            template<typename T>
            void
            HemmLLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Upper Hemm
            template<typename T>
            void
            HemmLU
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Upper Hemm where we avoid redistributing C
            template<typename T>
            void
            HemmLUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Lower Hemm
            template<typename T>
            void
            HemmRL
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Lower Hemm where we avoid redistributing C
            template<typename T>
            void
            HemmRLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Upper Hemm
            template<typename T>
            void
            HemmRU
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Upper Hemm where we avoid redistributing C
            template<typename T>
            void
            HemmRUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Her2k
            template<typename T>
            void
            Her2kLN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Her2k Update
            template<typename T>
            void
            Her2kLNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
                             const DistMatrix<T,MC,Star>& B_MC_Star,
                             const DistMatrix<T,MR,Star>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C       );

            // Lower, ConjugateTranspose Her2k
            template<typename T>
            void
            Her2kLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, ConjugateTranspose Her2k Update
            template<typename T>
            void
            Her2kLCUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_MC_Star,
                             const DistMatrix<T,Star,MR>& A_MR_Star,
                             const DistMatrix<T,Star,MC>& B_MC_Star,
                             const DistMatrix<T,Star,MR>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, Normal Her2k
            template<typename T>
            void
            Her2kUN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Normal Her2k Update
            template<typename T>
            void
            Her2kUNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
                             const DistMatrix<T,MC,Star>& B_MC_Star,
                             const DistMatrix<T,MR,Star>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C       );

            // Upper, ConjugateTranspose Her2k
            template<typename T>
            void
            Her2kUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, ConjugateTranspose Her2k Update
            template<typename T>
            void
            Her2kUCUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_MC_Star,
                             const DistMatrix<T,Star,MR>& A_MR_Star,
                             const DistMatrix<T,Star,MC>& B_MC_Star,
                             const DistMatrix<T,Star,MR>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,  MR>& C         );


            // Lower, Normal Herk
            template<typename T>
            void
            HerkLN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Herk Update
            template<typename T>
            void
            HerkLNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Lower, Normal Herk Update Kernel
            template<typename T>
            void
            HerkLNUpdateKernel
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Lower, ConjugateTranspose Herk
            template<typename T>
            void
            HerkLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, ConjugateTranspose Herk Update
            template<typename T>
            void
            HerkLCUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Lower, ConjugateTranspose Herk Update Kernel
            template<typename T>
            void
            HerkLCUpdateKernel
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, Normal Herk
            template<typename T>
            void
            HerkUN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Normal Herk Update
            template<typename T>
            void
            HerkUNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Upper, Normal Herk Update Kernel
            template<typename T>
            void
            HerkUNUpdateKernel
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Upper, ConjugateTranspose Herk
            template<typename T>
            void
            HerkUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, ConjugateTranspose Herk Update
            template<typename T>
            void
            HerkUCUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, ConjugateTranspose Herk Update Kernel
            template<typename T>
            void
            HerkUCUpdateKernel
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Symm
            template<typename T>
            void
            Symm
            ( const Side side, const Shape shape,
              const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Lower Symm
            template<typename T>
            void
            SymmLL
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Lower Symm where we avoid redistributing C
            template<typename T>
            void
            SymmLLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Upper Symm
            template<typename T>
            void
            SymmLU
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Left Upper Symm where we avoid redistributing C
            template<typename T>
            void
            SymmLUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Lower Symm
            template<typename T>
            void
            SymmRL
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Lower Symm where we avoid redistributing C
            template<typename T>
            void
            SymmRLC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Upper Symm
            template<typename T>
            void
            SymmRU
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Right Upper Symm where we avoid redistributing C
            template<typename T>
            void
            SymmRUC
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Syr2k
            template<typename T>
            void
            Syr2kLN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Syr2k Update
            template<typename T>
            void
            Syr2kLNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
                             const DistMatrix<T,MC,Star>& B_MC_Star,
                             const DistMatrix<T,MR,Star>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Lower, Transpose Syr2k
            template<typename T>
            void
            Syr2kLT
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Transpose Syr2k Update
            template<typename T>
            void
            Syr2kLTUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_MC_Star,
                             const DistMatrix<T,Star,MR>& A_MR_Star,
                             const DistMatrix<T,Star,MC>& B_MC_Star,
                             const DistMatrix<T,Star,MR>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, Normal Syr2k
            template<typename T>
            void
            Syr2kUN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Normal Syr2k Update
            template<typename T>
            void
            Syr2kUNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
                             const DistMatrix<T,MC,Star>& B_MC_Star,
                             const DistMatrix<T,MR,Star>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Upper, Transpose Syr2k
            template<typename T>
            void
            Syr2kUT
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
                             const DistMatrix<T,MC,MR>& B,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Transpose Syr2k Update
            template<typename T>
            void
            Syr2kUTUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_MC_Star,
                             const DistMatrix<T,Star,MR>& A_MR_Star,
                             const DistMatrix<T,Star,MC>& B_MC_Star,
                             const DistMatrix<T,Star,MR>& B_MR_Star,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Lower, Normal Syrk
            template<typename T>
            void
            SyrkLN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Normal Syrk Update
            template<typename T>
            void
            SyrkLNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Lower, Normal Syrk Update Kernel
            template<typename T>
            void
            SyrkLNUpdateKernel
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Lower, Transpose Syrk
            template<typename T>
            void
            SyrkLT
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Lower, Transpose Syrk Update
            template<typename T>
            void
            SyrkLTUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Lower, Transpose Syrk Update Kernel
            template<typename T>
            void
            SyrkLTUpdateKernel
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, Normal Syrk
            template<typename T>
            void
            SyrkUN
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Normal Syrk Update
            template<typename T>
            void
            SyrkUNUpdate
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Upper, Normal Syrk Update Kernel
            template<typename T>
            void
            SyrkUNUpdateKernel
            ( const T alpha, const DistMatrix<T,MC,Star>& A_MC_Star,
                             const DistMatrix<T,MR,Star>& A_MR_Star,
              const T beta,        DistMatrix<T,MC,MR  >& C         );

            // Upper, Transpose Syrk
            template<typename T>
            void
            SyrkUT
            ( const T alpha, const DistMatrix<T,MC,MR>& A,
              const T beta,        DistMatrix<T,MC,MR>& C );

            // Upper, Transpose Syrk Update
            template<typename T>
            void
            SyrkUTUpdate
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Upper, Transpose Syrk Update Kernel
            template<typename T>
            void
            SyrkUTUpdateKernel
            ( const T alpha, const DistMatrix<T,Star,MC>& A_Star_MC,
                             const DistMatrix<T,Star,MR>& A_Star_MR,
              const T beta,        DistMatrix<T,MC,  MR>& C         );

            // Left, Lower, Normal Trmm
            template<typename T>
            void
            TrmmLLN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Lower, (Conjugate)Transpose Trmm
            template<typename T>
            void
            TrmmLLT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Upper, Normal Trmm
            template<typename T>
            void
            TrmmLUN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Upper, (Conjugate)Transpose Trmm
            template<typename T>
            void
            TrmmLUT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Lower, Normal Trmm
            template<typename T>
            void
            TrmmRLN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Lower, (Conjugate)Transpose Trmm
            template<typename T>
            void
            TrmmRLT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Upper, Normal Trmm
            template<typename T>
            void
            TrmmRUN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Upper, (Conjugate)Transpose Trmm
            template<typename T>
            void
            TrmmRUT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Lower, Normal Trsm
            template<typename T>
            void
            TrsmLLN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Lower, (Conjugate)Transpose Trsm
            template<typename T>
            void
            TrsmLLT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Upper, Normal Trsm
            template<typename T>
            void
            TrsmLUN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Left, Upper, (Conjugate)Transpose Trsm
            template<typename T>
            void
            TrsmLUT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Lower, Normal Trsm
            template<typename T>
            void
            TrsmRLN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Lower, (Conjugate)Transpose Trsm
            template<typename T>
            void
            TrsmRLT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& L,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Upper, Normal Trsm
            template<typename T>
            void
            TrsmRUN
            ( const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            // Right, Upper, (Conjugate)Transpose Trsm
            template<typename T>
            void
            TrsmRUT
            ( const Orientation orientation, 
              const Diagonal diagonal,
              const T alpha, const DistMatrix<T,MC,MR>& U,
                                   DistMatrix<T,MC,MR>& X );

            //----------------------------------------------------------------//
            // Level 2 BLAS Utility Functions                                 //
            //----------------------------------------------------------------//
            template<typename T>
            double
            SymvGFlops
            ( const int m, const double seconds );
 
            //----------------------------------------------------------------//
            // Level 3 BLAS Utility Functions                                 //
            //----------------------------------------------------------------//
            template<typename T>
            double 
            GemmGFlops
            ( const int m, const int n, const int k, const double seconds );

            template<typename T>
            double
            HemmGFlops
            ( const Side side, const int m, const int n, const double seconds );

            template<typename T>
            double
            Her2kGFlops
            ( const int m, const int k, const double seconds );

            template<typename T>
            double
            HerkGFlops
            ( const int m, const int k, const double seconds );
            
            template<typename T>
            double
            SymmGFlops
            ( const Side side, const int m, const int n, const double seconds );
            
            template<typename T>
            double
            Syr2kGFlops
            ( const int m, const int k, const double seconds );

            template<typename T>
            double
            SyrkGFlops
            ( const int m, const int k, const double seconds );
 
            template<typename T>
            double
            TrmmGFlops
            ( const Side side, const int m, const int n, const double seconds );
            
            template<typename T>
            double
            TrsmGFlops
            ( const Side side, const int m, const int n, const double seconds );
        }
    }
}

/*----------------------------------------------------------------------------*/

namespace Elemental
{
    namespace BLAS
    {
        namespace Internal
        {

            // Level 2 Utility functions
            template<>
            inline double
            SymvGFlops<float>
            ( const int m, const double seconds )
            { return (1.*m*m)/(1.e9*seconds); }

            template<>
            inline double
            SymvGFlops<double>
            ( const int m, const double seconds )
            { return SymvGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            SymvGFlops<scomplex>
            ( const int m, const double seconds )
            { return 4.*SymvGFlops<float>(m,seconds); }

            template<>
            inline double
            SymvGFlops<dcomplex>
            ( const int m, const double seconds )
            { return 4.*SymvGFlops<float>(m,seconds); }
#endif

            // Level 3 Utility functions
            template<>
            inline double
            GemmGFlops<float>
            ( const int m, const int n, const int k, const double seconds )
            { return (2.*m*n*k)/(1.e9*seconds); }

            template<>
            inline double
            GemmGFlops<double>
            ( const int m, const int n, const int k, const double seconds )
            { return GemmGFlops<float>(m,n,k,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            GemmGFlops<scomplex>
            ( const int m, const int n, const int k, const double seconds )
            { return 4.*GemmGFlops<float>(m,n,k,seconds); }

            template<>
            inline double
            GemmGFlops<dcomplex>
            ( const int m, const int n, const int k, const double seconds )
            { return 4.*GemmGFlops<float>(m,n,k,seconds); }
#endif

            template<>
            inline double
            HemmGFlops<float>
            ( const Side side, const int m, const int n, const double seconds )
            {
                if( side == Left )
                    return (2.*m*m*n)/(1.e9*seconds);
                else
                    return (2.*m*n*n)/(1.e9*seconds);
            }

            template<>
            inline double
            HemmGFlops<double>
            ( const Side side, const int m, const int n, const double seconds )
            { return HemmGFlops<float>(side,m,n,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            HemmGFlops<scomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*HemmGFlops<float>(side,m,n,seconds); }

            template<>
            inline double
            HemmGFlops<dcomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*HemmGFlops<float>(side,m,n,seconds); }
#endif

            template<>
            inline double
            Her2kGFlops<float>
            ( const int m, const int k, const double seconds )
            { return (2.*m*m*k)/(1.e9*seconds); }

            template<>
            inline double
            Her2kGFlops<double>
            ( const int m, const int k, const double seconds )
            { return Her2kGFlops<float>(m,k,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            Her2kGFlops<scomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*Her2kGFlops<float>(m,k,seconds); }

            template<>
            inline double
            Her2kGFlops<dcomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*Her2kGFlops<float>(m,k,seconds); }
#endif

            template<>
            inline double
            HerkGFlops<float>
            ( const int m, const int k, const double seconds )
            { return (1.*m*m*k)/(1.e9*seconds); }

            template<>
            inline double
            HerkGFlops<double>
            ( const int m, const int k, const double seconds )
            { return HerkGFlops<float>(m,k,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            HerkGFlops<scomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*HerkGFlops<float>(m,k,seconds); }
            
            template<>
            inline double
            HerkGFlops<dcomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*HerkGFlops<float>(m,k,seconds); }
#endif
            
            template<>
            inline double
            SymmGFlops<float>
            ( const Side side, const int m, const int n, const double seconds )
            {
                if( side == Left )
                    return (2.*m*m*n)/(1.e9*seconds);
                else
                    return (2.*m*n*n)/(1.e9*seconds);
            }
            
            template<>
            inline double
            SymmGFlops<double>
            ( const Side side, const int m, const int n, const double seconds ) 
            { return SymmGFlops<float>(side,m,n,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            SymmGFlops<scomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*SymmGFlops<float>(side,m,n,seconds); }
            
            template<>
            inline double
            SymmGFlops<dcomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*SymmGFlops<float>(side,m,n,seconds); }
#endif
            
            template<>
            inline double
            Syr2kGFlops<float>
            ( const int m, const int k, const double seconds )
            { return (2.*m*m*k)/(1.e9*seconds); }
            
            template<>
            inline double
            Syr2kGFlops<double>
            ( const int m, const int k, const double seconds )
            { return Syr2kGFlops<float>(m,k,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            Syr2kGFlops<scomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*Syr2kGFlops<float>(m,k,seconds); }
            
            template<>
            inline double
            Syr2kGFlops<dcomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*Syr2kGFlops<float>(m,k,seconds); }
#endif
            
            template<>
            inline double
            SyrkGFlops<float>
            ( const int m, const int k, const double seconds )
            { return (1.*m*m*k)/(1.e9*seconds); }
            
            template<>
            inline double
            SyrkGFlops<double>
            ( const int m, const int k, const double seconds )
            { return SyrkGFlops<float>(m,k,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            SyrkGFlops<scomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*SyrkGFlops<float>(m,k,seconds); }
            
            template<>
            inline double
            SyrkGFlops<dcomplex>
            ( const int m, const int k, const double seconds )
            { return 4.*SyrkGFlops<float>(m,k,seconds); }
#endif
            
            template<>
            inline double
            TrmmGFlops<float>
            ( const Side side, const int m, const int n, const double seconds )
            {
                if( side == Left )
                    return (1.*m*m*n)/(1.e9*seconds);
                else
                    return (1.*m*n*n)/(1.e9*seconds);
            }
            
            template<>
            inline double
            TrmmGFlops<double>
            ( const Side side, const int m, const int n, const double seconds )
            { return TrmmGFlops<float>(side,m,n,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            TrmmGFlops<scomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*TrmmGFlops<float>(side,m,n,seconds); }
            
            template<>
            inline double
            TrmmGFlops<dcomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*TrmmGFlops<float>(side,m,n,seconds); }
#endif
            
            template<>
            inline double
            TrsmGFlops<float>
            ( const Side side, const int m, const int n, const double seconds )
            {
                if( side == Left )
                    return (1.*m*m*n)/(1.e9*seconds);
                else
                    return (1.*m*n*n)/(1.e9*seconds);
            }
            
            template<>
            inline double
            TrsmGFlops<double>
            ( const Side side, const int m, const int n, const double seconds )
            { return TrsmGFlops<float>(side,m,n,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            TrsmGFlops<scomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*TrsmGFlops<float>(side,m,n,seconds); }
            
            template<>
            inline double
            TrsmGFlops<dcomplex>
            ( const Side side, const int m, const int n, const double seconds )
            { return 4.*TrsmGFlops<float>(side,m,n,seconds); }
#endif
            
        }
    }
}

#endif /* ELEMENTAL_BLAS_INTERNAL_H */

