/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"
using namespace El;

#define RCDDM_s(AHandle) \
  reinterpret_cast<DynamicDistMatrix<float>*>(AHandle)
#define RCDDM_d(AHandle) \
  reinterpret_cast<DynamicDistMatrix<double>*>(AHandle)
#define RCDDM_c(AHandle) \
  reinterpret_cast<DynamicDistMatrix<Complex<float>>*>(AHandle)
#define RCDDM_z(AHandle) \
  reinterpret_cast<DynamicDistMatrix<Complex<double>>*>(AHandle)

#define RCDDM_s_const(AHandle) \
  reinterpret_cast<const DynamicDistMatrix<float>*>(AHandle)
#define RCDDM_d_const(AHandle) \
  reinterpret_cast<const DynamicDistMatrix<double>*>(AHandle)
#define RCDDM_c_const(AHandle) \
  reinterpret_cast<const DynamicDistMatrix<Complex<float>>*>(AHandle)
#define RCDDM_z_const(AHandle) \
  reinterpret_cast<const DynamicDistMatrix<Complex<double>>*>(AHandle)

#define RCADM_s(AHandle) \
  reinterpret_cast<AbstractDistMatrix<float>*>(AHandle)
#define RCADM_d(AHandle) \
  reinterpret_cast<AbstractDistMatrix<double>*>(AHandle)
#define RCADM_c(AHandle) \
  reinterpret_cast<AbstractDistMatrix<Complex<float>>*>(AHandle)
#define RCADM_z(AHandle) \
  reinterpret_cast<AbstractDistMatrix<Complex<double>>*>(AHandle)

#define RCADM_s_const(AHandle) \
  reinterpret_cast<AbstractDistMatrix<float>*>(AHandle)
#define RCADM_d_const(AHandle) \
  reinterpret_cast<AbstractDistMatrix<double>*>(AHandle)
#define RCADM_c_const(AHandle) \
  reinterpret_cast<AbstractDistMatrix<Complex<float>>*>(AHandle)
#define RCADM_z_const(AHandle) \
  reinterpret_cast<AbstractDistMatrix<Complex<double>>*>(AHandle)

#define RCDM_s(AHandle,ColDist,RowDist) \
  reinterpret_cast<DistMatrix<float,ColDist,RowDist>*>(AHandle)
#define RCDM_d(AHandle,ColDist,RowDist) \
  reinterpret_cast<DistMatrix<double,ColDist,RowDist>*>(AHandle)
#define RCDM_c(AHandle,ColDist,RowDist) \
  reinterpret_cast<DistMatrix<Complex<float>,ColDist,RowDist>*>(AHandle)
#define RCDM_z(AHandle,ColDist,RowDist) \
  reinterpret_cast<DistMatrix<Complex<double>,ColDist,RowDist>*>(AHandle)

#define RCDM_s_const(AHandle,ColDist,RowDist) \
  reinterpret_cast<const DistMatrix<float,ColDist,RowDist>*>(AHandle)
#define RCDM_d_const(AHandle,ColDist,RowDist) \
  reinterpret_cast<const DistMatrix<double,ColDist,RowDist>*>(AHandle)
#define RCDM_c_const(AHandle,ColDist,RowDist) \
  reinterpret_cast<const DistMatrix<Complex<float>,ColDist,RowDist>*>(AHandle)
#define RCDM_z_const(AHandle,ColDist,RowDist) \
  reinterpret_cast<const DistMatrix<Complex<double>,ColDist,RowDist>*>(AHandle)

#define DCDM_s(AHandle,ColDist,RowDist) \
  dynamic_cast<DistMatrix<float,ColDist,RowDist>*>(AHandle)
#define DCDM_d(AHandle,ColDist,RowDist) \
  dynamic_cast<DistMatrix<double,ColDist,RowDist>*>(AHandle)
#define DCDM_c(AHandle,ColDist,RowDist) \
  dynamic_cast<DistMatrix<Complex<float>,ColDist,RowDist>*>(AHandle)
#define DCDM_z(AHandle,ColDist,RowDist) \
  dynamic_cast<DistMatrix<Complex<double>,ColDist,RowDist>*>(AHandle)

#define DCDM_s_const(AHandle,ColDist,RowDist) \
  dynamic_cast<const DistMatrix<float,ColDist,RowDist>*>(AHandle)
#define DCDM_d_const(AHandle,ColDist,RowDist) \
  dynamic_cast<const DistMatrix<double,ColDist,RowDist>*>(AHandle)
#define DCDM_c_const(AHandle,ColDist,RowDist) \
  dynamic_cast<const DistMatrix<Complex<float>,ColDist,RowDist>*>(AHandle)
#define DCDM_z_const(AHandle,ColDist,RowDist) \
  dynamic_cast<const DistMatrix<Complex<double>,ColDist,RowDist>*>(AHandle)

#define RCB_c(buffer) reinterpret_cast<Complex<float>*>(buffer)
#define RCB_z(buffer) reinterpret_cast<Complex<double>*>(buffer)

#define RCB_c_const(buffer) reinterpret_cast<const Complex<float>*>(buffer)
#define RCB_z_const(buffer) reinterpret_cast<const Complex<double>*>(buffer)

#define DC_CHECK(A) if( A == nullptr ) RuntimeError("Dynamic cast failed");

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// DistMatrix::Get( Int i, Int j ) const
// -------------------------------------

float ElDistMatrixGet_s( const ElDistMatrix_s* AHandle, ElInt i, ElInt j )
{
    float alpha;
    try { alpha = RCDDM_s_const(AHandle)->ADM->Get(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGet_d( const ElDistMatrix_d* AHandle, ElInt i, ElInt j )
{
    double alpha;
    try { alpha = RCDDM_d_const(AHandle)->ADM->Get(i,j); }
    CATCH
    return alpha;
}

void ElDistMatrixGet_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_c(alpha) = RCDDM_c_const(AHandle)->ADM->Get(i,j); }
    CATCH
}

void ElDistMatrixGet_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j, void* alpha )
{
    try { *RCB_z(alpha) = RCDDM_z_const(AHandle)->ADM->Get(i,j); }
    CATCH
}

// B = A
// -----

#define INNER_IF_CONVERT_AND_COPY(A,BDyn,ColDist,RowDist) \
  if( BDyn->U == ColDist && BDyn->V == RowDist ) \
  { \
      auto B = DCDM_s(BDyn->ADM,ColDist,RowDist); \
      DC_CHECK(B); \
      *B = *A; \
  }
#define INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,ColDist,RowDist) \
  else INNER_IF_CONVERT_AND_COPY(A,BDyn,ColDist,RowDist)

#define IF_CONVERT_AND_COPY(ADyn,BDyn,ColDist,RowDist) \
  if( ADyn->U == ColDist && ADyn->V == RowDist ) \
  { \
      auto A = DCDM_s_const(ADyn->ADM,ColDist,RowDist); \
      DC_CHECK(A); \
      INNER_IF_CONVERT_AND_COPY(A,BDyn,CIRC,CIRC) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,MC,  MR  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,MC,  STAR) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,MD,  STAR) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,MR,  MC  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,MR,  STAR) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,MC  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,MD  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,MR  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,STAR) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,VC  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,STAR,VR  ) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,VC,  STAR) \
      INNER_ELSEIF_CONVERT_AND_COPY(A,BDyn,VR,  STAR) \
  }
#define ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,ColDist,RowDist) \
  else IF_CONVERT_AND_COPY(ADyn,BDyn,ColDist,RowDist)

void ElDistMatrixCopy_s
( const ElDistMatrix_s* AHandle, ElDistMatrix_s* BHandle )
{
    try
    {
        auto ADyn = RCDDM_s_const(AHandle);
        auto BDyn = RCDDM_s(BHandle);
        IF_CONVERT_AND_COPY(ADyn,BDyn,CIRC,CIRC)
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MC,  MR  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MC,  STAR) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MD,  STAR) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MR,  MC  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MR,  STAR) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,MC  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,MD  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,MR  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,STAR) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,VC  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,VR  ) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,VC,  STAR) 
        ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,VR,  STAR) 
    }
    CATCH
}

} // extern "C"
