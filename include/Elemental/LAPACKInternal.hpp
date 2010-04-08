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
#ifndef ELEMENTAL_LAPACK_INTERNAL_HPP
#define ELEMENTAL_LAPACK_INTERNAL_HPP 1

#include "Elemental/LAPACK.hpp"

namespace Elemental
{
    namespace LAPACK
    {
        namespace Internal
        {
            //----------------------------------------------------------------//
            // Chol                                                           //
            //----------------------------------------------------------------//
            template<typename T>
            void
            CholVar2
            ( const Shape shape, DistMatrix<T,MC,MR>& A );

            template<typename T>
            void
            CholVar3
            ( const Shape shape, DistMatrix<T,MC,MR>& A );

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
            
            //----------------------------------------------------------------//
            // GaussElim                                                      //
            //----------------------------------------------------------------//
            
            template<typename T>
            void
            ReduceToRowEchelon
            ( DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B );

            //----------------------------------------------------------------//
            // LU                                                             //
            //----------------------------------------------------------------//

            template<typename T>
            void
            ApplyRowPivots
            (       DistMatrix<T,MC,MR>& A, 
              const std::vector<int>& image,
              const std::vector<int>& preimage,
              const int pivotOffset=0           );

            void
            ComposePivots
            ( const DistMatrix<int,Star,Star>& p,
                    std::vector<int>& image,
                    std::vector<int>& preimage,
              const int pivotOffset = 0          );

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
              const int pivotOffset=0      );

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

            //----------------------------------------------------------------//
            // Tridiag                                                        //
            //----------------------------------------------------------------//
            template<typename T>
            T
            Reflector( DistMatrix<T,MC,MR>& x );

            template<typename T>
            T
            LocalColReflector( DistMatrix<T,MC,MR>& x );

            template<typename T>
            T
            LocalRowReflector( DistMatrix<T,MC,MR>& x );

            template<typename T>
            void
            PanelTridiagL
            ( DistMatrix<T,MC,  MR  >& A,
              DistMatrix<T,MC,  MR  >& W,
              DistMatrix<T,MD,  Star>& e,
              DistMatrix<T,MD,  Star>& t,
              DistMatrix<T,MC,  Star>& A_MC_Star,
              DistMatrix<T,Star,MR  >& ATrans_Star_MR );
 
            template<typename T>
            void
            PanelTridiagU
            ( DistMatrix<T,MC,MR  >& A,
              DistMatrix<T,MC,MR  >& W,
              DistMatrix<T,MD,Star>& e,
              DistMatrix<T,MD,Star>& t );

            template<typename T>
            void
            TridiagL
            ( DistMatrix<T,MC,MR  >& A,
              DistMatrix<T,MD,Star>& d,
              DistMatrix<T,MD,Star>& e,
              DistMatrix<T,MD,Star>& t );
           
            template<typename T>
            void
            TridiagU
            ( DistMatrix<T,MC,MR  >& A,
              DistMatrix<T,MD,Star>& d,
              DistMatrix<T,MD,Star>& e,
              DistMatrix<T,MD,Star>& t );

            //----------------------------------------------------------------//
            // Trinv                                                          //
            //----------------------------------------------------------------//
            template<typename T>
            void
            TrinvVar3
            ( const Shape shape, 
              const Diagonal diagonal,
              DistMatrix<T,MC,MR>& A  );

            template<typename T>
            void
            TrinvL
            ( const Diagonal diagonal, DistMatrix<T,MC,MR>& L );

            template<typename T>
            void
            TrinvU
            ( const Diagonal diagonal, DistMatrix<T,MC,MR>& U );

            template<typename T>
            void
            TrinvLVar3
            ( const Diagonal diagonal, DistMatrix<T,MC,MR>& L );

            template<typename T>
            void
            TrinvUVar3
            ( const Diagonal diagonal, DistMatrix<T,MC,MR>& U );

            //----------------------------------------------------------------//
            // LAPACK Utility Functions                                       //
            //----------------------------------------------------------------//
            template<typename T>
            double
            CholGFlops( const int m, const double seconds );

            template<typename T>
            double
            LUGFlops( const int m, const double seconds );

            template<typename T>
            double
            TridiagGFlops( const int m, const double seconds );

            template<typename T>
            double
            TrinvGFlops( const int m, const double seconds );
        }
    }
}

/*----------------------------------------------------------------------------*/

namespace Elemental
{
    namespace LAPACK
    {
        namespace Internal
        {

            template<>
            inline double
            CholGFlops<float>
            ( const int m, const double seconds )
            { return (1./3.*m*m*m)/(1.e9*seconds); }
            
            template<>
            inline double
            CholGFlops<double>
            ( const int m, const double seconds )
            { return CholGFlops<float>(m,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            CholGFlops<scomplex>
            ( const int m, const double seconds )
            { return 4.*CholGFlops<float>(m,seconds); }
            
            template<>
            inline double
            CholGFlops<dcomplex>
            ( const int m, const double seconds )
            { return 4.*CholGFlops<float>(m,seconds); }
#endif

            template<>
            inline double
            LUGFlops<float>
            ( const int m, const double seconds )
            { return (2./3.*m*m*m)/(1.e9*seconds); }

            template<>
            inline double
            LUGFlops<double>
            ( const int m, const double seconds )
            { return LUGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            LUGFlops<scomplex>
            ( const int m, const double seconds )
            { return 4.*LUGFlops<float>(m,seconds); }

            template<>
            inline double
            LUGFlops<dcomplex>
            ( const int m, const double seconds )
            { return 4.*LUGFlops<float>(m,seconds); }
#endif

            template<>
            inline double
            TridiagGFlops<float>
            ( const int m, const double seconds )
            { return (4./3.*m*m*m)/(1.e9*seconds); }

            template<>
            inline double
            TridiagGFlops<double>
            ( const int m, const double seconds )
            { return TridiagGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            TridiagGFlops<scomplex>
            ( const int m, const double seconds )
            { return 4.*TridiagGFlops<float>(m,seconds); }

            template<>
            inline double
            TridiagGFlops<dcomplex>
            ( const int m, const double seconds )
            { return 4.*TridiagGFlops<float>(m,seconds); }
#endif

            template<>
            inline double
            TrinvGFlops<float>
            ( const int m, const double seconds )
            { return (1./3.*m*m*m)/(1.e9*seconds); }
            
            template<>
            inline double
            TrinvGFlops<double>
            ( const int m, const double seconds )
            { return TrinvGFlops<float>(m,seconds); }
            
#ifndef WITHOUT_COMPLEX
            template<>
            inline double
            TrinvGFlops<scomplex>
            ( const int m, const double seconds )
            { return 4.*TrinvGFlops<float>(m,seconds); }
            
            template<>
            inline double
            TrinvGFlops<dcomplex>
            ( const int m, const double seconds )
            { return 4.*TrinvGFlops<float>(m,seconds); }
#endif

        }
    }
}

#endif /* ELEMENTAL_LAPACK_INTERNAL_HPP */

