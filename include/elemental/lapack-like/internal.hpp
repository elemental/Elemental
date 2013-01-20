/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_INTERNAL_HPP
#define LAPACK_INTERNAL_HPP

namespace elem {
namespace internal {

//----------------------------------------------------------------------------//
// Local LAPACK-like routines                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalLDL
( Orientation orientation, 
  DistMatrix<F,STAR,STAR>& A, DistMatrix<F,STAR,STAR>& d );

template<typename F>
void LocalLU( DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );

//----------------------------------------------------------------------------//
// Cholesky helpers                                                           //
//----------------------------------------------------------------------------//

template<typename F>
void CholeskyLVar3Unb( Matrix<F>& A );
template<typename F>
void CholeskyLVar2( Matrix<F>& A );
template<typename F>
void CholeskyLVar2( DistMatrix<F>& A );
template<typename F>
void CholeskyLVar3( Matrix<F>& A );
template<typename F>
void CholeskyLVar3( DistMatrix<F>& A );
template<typename F>
void CholeskyLVar3Square( DistMatrix<F>& A );

template<typename F>
void CholeskyUVar3Unb( Matrix<F>& A );
template<typename F>
void CholeskyUVar2( Matrix<F>& A );
template<typename F>
void CholeskyUVar2( DistMatrix<F>& A );
template<typename F>
void CholeskyUVar3( Matrix<F>& A );
template<typename F>
void CholeskyUVar3( DistMatrix<F>& A );
template<typename F>
void CholeskyUVar3Square( DistMatrix<F>& A );
            
//----------------------------------------------------------------------------//
// GaussElim                                                                  //
//----------------------------------------------------------------------------//
            
template<typename F>
void ReduceToRowEchelon( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void ReduceToRowEchelon( DistMatrix<F>& A, DistMatrix<F>& B );

//----------------------------------------------------------------------------//
// LDL                                                                        //
//----------------------------------------------------------------------------//

template<typename F>
void LDLVar3Unb
( Orientation orientation, Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLVar3
( Orientation orientation, Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLVar3
( Orientation orientation, DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d );

//----------------------------------------------------------------------------//
// LU                                                                         //
//----------------------------------------------------------------------------//

void ComposePanelPivots
( const Matrix<int>& p,
        int pivotOffset,
        std::vector<int>& image,
        std::vector<int>& preimage );

void ComposePanelPivots
( const DistMatrix<int,STAR,STAR>& p,
        int pivotOffset,
        std::vector<int>& image,
        std::vector<int>& preimage );

template<typename F> void CreatePivotOp();
template<> void CreatePivotOp<float>();
template<> void CreatePivotOp<double>();
template<> void CreatePivotOp<scomplex>();
template<> void CreatePivotOp<dcomplex>();

template<typename T> void DestroyPivotOp();
template<> void DestroyPivotOp<float>();
template<> void DestroyPivotOp<double>();
template<> void DestroyPivotOp<scomplex>();
template<> void DestroyPivotOp<dcomplex>();

template<typename F>
void PanelLU( Matrix<F>& A, Matrix<int>& p, int pivotOffset=0 );
template<typename F>
void PanelLU
( DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<F,MC,  STAR>& B, 
  DistMatrix<int,STAR,STAR>& p, 
  int pivotOffset=0 );

template<typename F> mpi::Op PivotOp();
template<> mpi::Op PivotOp<float>();
template<> mpi::Op PivotOp<double>();
template<> mpi::Op PivotOp<scomplex>();
template<> mpi::Op PivotOp<dcomplex>();
            
template<typename F>
void PivotFunc
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

//----------------------------------------------------------------------------//
// LQ                                                                         //
//----------------------------------------------------------------------------//

template<typename R>
void PanelLQ( Matrix<R>& A );
template<typename R>
void PanelLQ( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );

template<typename R>
void PanelLQ( DistMatrix<R>& A );
template<typename R>
void PanelLQ
( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t );

//----------------------------------------------------------------------------//
// QR                                                                         //
//----------------------------------------------------------------------------//

template<typename R>
void PanelQR( Matrix<R>& A );
template<typename R>
void PanelQR( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );

template<typename R>
void PanelQR( DistMatrix<R>& A );
template<typename R>
void PanelQR
( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t );

//----------------------------------------------------------------------------//
// Reflector                                                                  //
//----------------------------------------------------------------------------//
template<typename R>
R ColReflector( DistMatrix<R>& chi, DistMatrix<R>& x );

template<typename R>
Complex<R> ColReflector
( DistMatrix<Complex<R> >& chi, 
  DistMatrix<Complex<R> >& x );

template<typename R>
R RowReflector( DistMatrix<R>& chi, DistMatrix<R>& x );

template<typename R>
Complex<R> RowReflector
( DistMatrix<Complex<R> >& chi,
  DistMatrix<Complex<R> >& x );

//----------------------------------------------------------------------------//
// Bidiag                                                                     //
//----------------------------------------------------------------------------//

template<typename R>
void BidiagL( DistMatrix<R>& A );
template<typename R>
void BidiagU( DistMatrix<R>& A );

template<typename R>
void BidiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ );
template<typename R>
void BidiagU
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ );

template<typename R>
void PanelBidiagL
( DistMatrix<R>& A,
  DistMatrix<R>& X,
  DistMatrix<R>& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR );
template<typename R>
void PanelBidiagU
( DistMatrix<R>& A,
  DistMatrix<R>& X,
  DistMatrix<R>& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR );

template<typename R>
void PanelBidiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,MD,  STAR>& tP,
  DistMatrix<Complex<R>,MD,  STAR>& tQ,
  DistMatrix<Complex<R> >& X,
  DistMatrix<Complex<R> >& Y,
  DistMatrix<Complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<Complex<R>,STAR,MR  >& ARowPan_STAR_MR );
template<typename R>
void PanelBidiagU
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,MD,  STAR>& tP,
  DistMatrix<Complex<R>,MD,  STAR>& tQ,
  DistMatrix<Complex<R> >& X,
  DistMatrix<Complex<R> >& Y,
  DistMatrix<Complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<Complex<R>,STAR,MR  >& ARowPan_STAR_MR );

//----------------------------------------------------------------------------//
// HermitianTridiag                                                           //
//----------------------------------------------------------------------------//

template<typename R>
void HermitianPanelTridiagL
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagU
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagLSquare
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagUSquare
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );

template<typename R>
void HermitianPanelTridiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagU
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );

template<typename R>
void HermitianPanelTridiagLSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagUSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
 
template<typename R>
void HermitianTridiagL( DistMatrix<R>& A );
template<typename R>
void HermitianTridiagU( DistMatrix<R>& A );

template<typename R>
void HermitianTridiagLSquare( DistMatrix<R>& A );
template<typename R>
void HermitianTridiagUSquare( DistMatrix<R>& A );

template<typename R>
void HermitianTridiagL
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );
template<typename R>
void HermitianTridiagU
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );

template<typename R>
void HermitianTridiagLSquare
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );
template<typename R>
void HermitianTridiagUSquare
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );

//----------------------------------------------------------------------------//
// HPD Inverse                                                                //
//----------------------------------------------------------------------------//

template<typename F>
void HPDInverseLVar2( DistMatrix<F>& A );
template<typename F>
void HPDInverseUVar2( DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// Triangular Inverse                                                         //
//----------------------------------------------------------------------------//

template<typename F>
void TriangularInverseVar3
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A  );
template<typename F>
void TriangularInverseVar3
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  );

template<typename F>
void TriangularInverseLVar3
( UnitOrNonUnit diag, Matrix<F>& L );
template<typename F>
void TriangularInverseLVar3
( UnitOrNonUnit diag, DistMatrix<F>& L );

template<typename F>
void TriangularInverseUVar3
( UnitOrNonUnit diag, Matrix<F>& U );
template<typename F>
void TriangularInverseUVar3
( UnitOrNonUnit diag, DistMatrix<F>& U );

} // namespace internal
} // namespace elem

#endif // ifndef LAPACK_INTERNAL_HPP
