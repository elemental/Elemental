/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.
   Copyright (c) 2014, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include EL_COPY_INC 

#include "El-C.h"
using namespace El;

#define RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

#define RCG(grid) RC(Grid*,grid)
#define RCG_const(grid) RC(const Grid*,grid)

#define RCADM_s(A) RC(AbstractDistMatrix<float          >*,A)
#define RCADM_d(A) RC(AbstractDistMatrix<double         >*,A)
#define RCADM_c(A) RC(AbstractDistMatrix<Complex<float >>*,A)
#define RCADM_z(A) RC(AbstractDistMatrix<Complex<double>>*,A)

#define RCADM_s_const(A) RC(const AbstractDistMatrix<float          >*,A)
#define RCADM_d_const(A) RC(const AbstractDistMatrix<double         >*,A)
#define RCADM_c_const(A) RC(const AbstractDistMatrix<Complex<float >>*,A)
#define RCADM_z_const(A) RC(const AbstractDistMatrix<Complex<double>>*,A)

#define RCB_c(buffer) RC(Complex<float>*,buffer)
#define RCB_z(buffer) RC(Complex<double>*,buffer)

#define RCB_c_const(buffer) RC(const Complex<float >*,buffer)
#define RCB_z_const(buffer) RC(const Complex<double>*,buffer)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

// DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
// ------------------------------------------------
ElDistMatrix_s* ElDistMatrixCreate_s( const ElGrid gridHandle )
{
    ElDistMatrix_s* AHandle = 0;
    try 
    {
        const Grid* grid = RCG_const(gridHandle);
        AHandle = RC(ElDistMatrix_s*,new DistMatrix<float>(*grid));
    }
    CATCH
    return AHandle;
}

ElDistMatrix_d* ElDistMatrixCreate_d( const ElGrid gridHandle )
{
    ElDistMatrix_d* AHandle = 0;
    try 
    {
        const Grid* grid = RCG_const(gridHandle);
        AHandle = RC(ElDistMatrix_d*,new DistMatrix<double>(*grid));
    }
    CATCH
    return AHandle;
}

ElDistMatrix_c* ElDistMatrixCreate_c( const ElGrid gridHandle )
{
    ElDistMatrix_c* AHandle = 0;
    try 
    {
        const Grid* grid = RCG_const(gridHandle);
        AHandle = RC(ElDistMatrix_c*,new DistMatrix<Complex<float>>(*grid));
    }
    CATCH
    return AHandle;
}

ElDistMatrix_z* ElDistMatrixCreate_z( const ElGrid gridHandle )
{
    ElDistMatrix_z* AHandle = 0;
    try 
    {
        const Grid* grid = RCG_const(gridHandle);
        AHandle = RC(ElDistMatrix_z*,new DistMatrix<Complex<double>>(*grid));
    }
    CATCH
    return AHandle;
}

// DistMatrix<T,U,V>::DistMatrix( const Grid& g )
// ----------------------------------------------
ElDistMatrix_s* ElDistMatrixCreateSpecific_s
( ElDist U_C, ElDist V_C, const ElGrid gridHandle )
{
    ElDistMatrix_s* AHandle = 0;
    try 
    {
        Dist U = static_cast<Dist>(U_C);
        Dist V = static_cast<Dist>(V_C);
        const Grid* grid = RCG_const(gridHandle);

        AbstractDistMatrix<float>* ADM;
        if( U == CIRC && V == CIRC )
            ADM = new DistMatrix<float,CIRC,CIRC>(*grid); 
        else if( U == MC && V == MR )
            ADM = new DistMatrix<float,MC,MR>(*grid);
        else if( U == MC && V == STAR )
            ADM = new DistMatrix<float,MC,STAR>(*grid);
        else if( U == MD && V == STAR )
            ADM = new DistMatrix<float,MD,STAR>(*grid);
        else if( U == MR && V == MC )
            ADM = new DistMatrix<float,MR,MC>(*grid);
        else if( U == MR && V == STAR )
            ADM = new DistMatrix<float,MR,STAR>(*grid);
        else if( U == STAR && V == MC )
            ADM = new DistMatrix<float,STAR,MC>(*grid);
        else if( U == STAR && V == MD )
            ADM = new DistMatrix<float,STAR,MD>(*grid);
        else if( U == STAR && V == MR )
            ADM = new DistMatrix<float,STAR,MR>(*grid);
        else if( U == STAR && V == STAR )
            ADM = new DistMatrix<float,STAR,STAR>(*grid);
        else if( U == STAR && V == VC )
            ADM = new DistMatrix<float,STAR,VC>(*grid);
        else if( U == STAR && V == VR )
            ADM = new DistMatrix<float,STAR,VR>(*grid);
        else if( U == VC && V == STAR )
            ADM = new DistMatrix<float,VC,STAR>(*grid);
        else if( U == VR && V == STAR )
            ADM = new DistMatrix<float,VR,STAR>(*grid);
        else
            RuntimeError("Invalid distribution pair");

        AHandle = RC(ElDistMatrix_s*,ADM);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_d* ElDistMatrixCreateSpecific_d
( ElDist U_C, ElDist V_C, const ElGrid gridHandle )
{
    ElDistMatrix_d* AHandle = 0;
    try 
    {
        Dist U = static_cast<Dist>(U_C);
        Dist V = static_cast<Dist>(V_C);
        const Grid* grid = RCG_const(gridHandle);

        AbstractDistMatrix<double>* ADM;
        if( U == CIRC && V == CIRC )
            ADM = new DistMatrix<double,CIRC,CIRC>(*grid);
        else if( U == MC && V == MR )
            ADM = new DistMatrix<double,MC,MR>(*grid);
        else if( U == MC && V == STAR )
            ADM = new DistMatrix<double,MC,STAR>(*grid);
        else if( U == MD && V == STAR )
            ADM = new DistMatrix<double,MD,STAR>(*grid);
        else if( U == MR && V == MC )
            ADM = new DistMatrix<double,MR,MC>(*grid);
        else if( U == MR && V == STAR )
            ADM = new DistMatrix<double,MR,STAR>(*grid);
        else if( U == STAR && V == MC )
            ADM = new DistMatrix<double,STAR,MC>(*grid);
        else if( U == STAR && V == MD )
            ADM = new DistMatrix<double,STAR,MD>(*grid);
        else if( U == STAR && V == MR )
            ADM = new DistMatrix<double,STAR,MR>(*grid);
        else if( U == STAR && V == STAR )
            ADM = new DistMatrix<double,STAR,STAR>(*grid);
        else if( U == STAR && V == VC )
            ADM = new DistMatrix<double,STAR,VC>(*grid);
        else if( U == STAR && V == VR )
            ADM = new DistMatrix<double,STAR,VR>(*grid);
        else if( U == VC && V == STAR )
            ADM = new DistMatrix<double,VC,STAR>(*grid);
        else if( U == VR && V == STAR )
            ADM = new DistMatrix<double,VR,STAR>(*grid);
        else
            RuntimeError("Invalid distribution pair");

        AHandle = RC(ElDistMatrix_d*,ADM);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_c* ElDistMatrixCreateSpecific_c
( ElDist U_C, ElDist V_C, const ElGrid gridHandle )
{
    ElDistMatrix_c* AHandle = 0;
    try 
    {
        Dist U = static_cast<Dist>(U_C);
        Dist V = static_cast<Dist>(V_C);
        const Grid* grid = RCG_const(gridHandle);

        AbstractDistMatrix<Complex<float>>* ADM;
        if( U == CIRC && V == CIRC )
            ADM = new DistMatrix<Complex<float>,CIRC,CIRC>(*grid);
        else if( U == MC && V == MR )
            ADM = new DistMatrix<Complex<float>,MC,MR>(*grid);
        else if( U == MC && V == STAR )
            ADM = new DistMatrix<Complex<float>,MC,STAR>(*grid);
        else if( U == MD && V == STAR )
            ADM = new DistMatrix<Complex<float>,MD,STAR>(*grid);
        else if( U == MR && V == MC )
            ADM = new DistMatrix<Complex<float>,MR,MC>(*grid);
        else if( U == MR && V == STAR )
            ADM = new DistMatrix<Complex<float>,MR,STAR>(*grid);
        else if( U == STAR && V == MC )
            ADM = new DistMatrix<Complex<float>,STAR,MC>(*grid);
        else if( U == STAR && V == MD )
            ADM = new DistMatrix<Complex<float>,STAR,MD>(*grid);
        else if( U == STAR && V == MR )
            ADM = new DistMatrix<Complex<float>,STAR,MR>(*grid);
        else if( U == STAR && V == STAR )
            ADM = new DistMatrix<Complex<float>,STAR,STAR>(*grid);
        else if( U == STAR && V == VC )
            ADM = new DistMatrix<Complex<float>,STAR,VC>(*grid);
        else if( U == STAR && V == VR )
            ADM = new DistMatrix<Complex<float>,STAR,VR>(*grid);
        else if( U == VC && V == STAR )
            ADM = new DistMatrix<Complex<float>,VC,STAR>(*grid);
        else if( U == VR && V == STAR )
            ADM = new DistMatrix<Complex<float>,VR,STAR>(*grid);
        else
            RuntimeError("Invalid distribution pair");

        AHandle = RC(ElDistMatrix_c*,ADM);
    }
    CATCH
    return AHandle;
}

ElDistMatrix_z* ElDistMatrixCreateSpecific_z
( ElDist U_C, ElDist V_C, const ElGrid gridHandle )
{
    ElDistMatrix_z* AHandle = 0;
    try 
    {
        Dist U = static_cast<Dist>(U_C);
        Dist V = static_cast<Dist>(V_C);
        const Grid* grid = RCG_const(gridHandle);

        AbstractDistMatrix<Complex<double>>* ADM;
        if( U == CIRC && V == CIRC )
            ADM = new DistMatrix<Complex<double>,CIRC,CIRC>(*grid);
        else if( U == MC && V == MR ) 
            ADM = new DistMatrix<Complex<double>,MC,MR>(*grid);
        else if( U == MC && V == STAR )
            ADM = new DistMatrix<Complex<double>,MC,STAR>(*grid);
        else if( U == MD && V == STAR )
            ADM = new DistMatrix<Complex<double>,MD,STAR>(*grid);
        else if( U == MR && V == MC )
            ADM = new DistMatrix<Complex<double>,MR,MC>(*grid);
        else if( U == MR && V == STAR )
            ADM = new DistMatrix<Complex<double>,MR,STAR>(*grid);
        else if( U == STAR && V == MC )
            ADM = new DistMatrix<Complex<double>,STAR,MC>(*grid);
        else if( U == STAR && V == MD )
            ADM = new DistMatrix<Complex<double>,STAR,MD>(*grid);
        else if( U == STAR && V == MR )
            ADM = new DistMatrix<Complex<double>,STAR,MR>(*grid);
        else if( U == STAR && V == STAR )
            ADM = new DistMatrix<Complex<double>,STAR,STAR>(*grid);
        else if( U == STAR && V == VC )
            ADM = new DistMatrix<Complex<double>,STAR,VC>(*grid);
        else if( U == STAR && V == VR )
            ADM = new DistMatrix<Complex<double>,STAR,VR>(*grid);
        else if( U == VC && V == STAR )
            ADM = new DistMatrix<Complex<double>,VC,STAR>(*grid);
        else if( U == VR && V == STAR )
            ADM = new DistMatrix<Complex<double>,VR,STAR>(*grid);
        else
            RuntimeError("Invalid distribution pair");

        AHandle = RC(ElDistMatrix_z*,ADM);
    }
    CATCH
    return AHandle;
}

// DistMatrix<T,U,V>::~DistMatrix()
// --------------------------------
void ElDistMatrixDestroy_s( const ElDistMatrix_s* AHandle )
{ delete RCADM_s_const(AHandle); }

void ElDistMatrixDestroy_d( const ElDistMatrix_d* AHandle )
{ delete RCADM_d_const(AHandle); }

void ElDistMatrixDestroy_c( const ElDistMatrix_c* AHandle )
{ delete RCADM_c_const(AHandle); }

void ElDistMatrixDestroy_z( const ElDistMatrix_z* AHandle )
{ delete RCADM_z_const(AHandle); }

// void DistMatrix<T,U,V>::Empty()
// -------------------------------
void ElDistMatrixEmpty_s( ElDistMatrix_s* AHandle )
{
    try { RCADM_s(AHandle)->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_d( ElDistMatrix_d* AHandle )
{
    try { RCADM_d(AHandle)->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_c( ElDistMatrix_c* AHandle )
{
    try { RCADM_c(AHandle)->Empty(); }
    CATCH
}

void ElDistMatrixEmpty_z( ElDistMatrix_z* AHandle )
{
    try { RCADM_z(AHandle)->Empty(); }
    CATCH
}

// void DistMatrix<T,U,V>::EmptyData()
// -----------------------------------
void ElDistMatrixEmptyData_s( ElDistMatrix_s* AHandle )
{
    try { RCADM_s(AHandle)->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_d( ElDistMatrix_d* AHandle )
{
    try { RCADM_d(AHandle)->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_c( ElDistMatrix_c* AHandle )
{
    try { RCADM_c(AHandle)->EmptyData(); }
    CATCH
}

void ElDistMatrixEmptyData_z( ElDistMatrix_z* AHandle )
{
    try { RCADM_z(AHandle)->EmptyData(); }
    CATCH
}

// void DistMatrix<T,U,V>::SetGrid( const Grid& g )
// ------------------------------------------------
void ElDistMatrixSetGrid_s( ElDistMatrix_s* AHandle, const ElGrid gridHandle )
{
    try { RCADM_s(AHandle)->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_d( ElDistMatrix_d* AHandle, const ElGrid gridHandle )
{
    try { RCADM_d(AHandle)->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_c( ElDistMatrix_c* AHandle, const ElGrid gridHandle )
{
    try { RCADM_c(AHandle)->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

void ElDistMatrixSetGrid_z( ElDistMatrix_z* AHandle, const ElGrid gridHandle )
{
    try { RCADM_z(AHandle)->SetGrid(*RCG_const(gridHandle)); }
    CATCH
}

// B = A
// -----
void ElDistMatrixCopy_s
( const ElDistMatrix_s* AHandle, ElDistMatrix_s* BHandle )
{ 
    try { Copy( *RCADM_s_const(AHandle), *RCADM_s(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_d
( const ElDistMatrix_d* AHandle, ElDistMatrix_d* BHandle )
{
    try { Copy( *RCADM_d_const(AHandle), *RCADM_d(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_c
( const ElDistMatrix_c* AHandle, ElDistMatrix_c* BHandle )
{
    try { Copy( *RCADM_c_const(AHandle), *RCADM_c(BHandle) ); }
    CATCH
}

void ElDistMatrixCopy_z
( const ElDistMatrix_z* AHandle, ElDistMatrix_z* BHandle )
{
    try { Copy( *RCADM_z_const(AHandle), *RCADM_z(BHandle) ); }
    CATCH
}

// void DistMatrix<T,U,V>::Resize( Int height, Int width )
// -------------------------------------------------------
void ElDistMatrixResize_s( ElDistMatrix_s* AHandle, ElInt height, ElInt width )
{
    try { RCADM_s(AHandle)->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_d( ElDistMatrix_d* AHandle, ElInt height, ElInt width )
{
    try { RCADM_d(AHandle)->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_c( ElDistMatrix_c* AHandle, ElInt height, ElInt width )
{
    try { RCADM_c(AHandle)->Resize(height,width); }
    CATCH
}

void ElDistMatrixResize_z( ElDistMatrix_z* AHandle, ElInt height, ElInt width )
{
    try { RCADM_z(AHandle)->Resize(height,width); }
    CATCH
}

// void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
// -----------------------------------------------------------------
void ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCADM_s(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCADM_d(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCADM_c(AHandle)->Resize(height,width,ldim); }
    CATCH
}

void ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z* AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { RCADM_z(AHandle)->Resize(height,width,ldim); }
    CATCH
}

// void DistMatrix<T,U,V>::MakeConsistent()
// ----------------------------------------
void ElDistMatrixMakeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers )
{
    try { RCADM_s(AHandle)->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers )
{
    try { RCADM_d(AHandle)->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers )
{
    try { RCADM_c(AHandle)->MakeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers )
{
    try { RCADM_z(AHandle)->MakeConsistent(includeViewers); }
    CATCH
}

// void DistMatrix<T,U,V>::MakeSizeConsistent()
// --------------------------------------------
void ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s* AHandle, bool includeViewers )
{
    try { RCADM_s(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d* AHandle, bool includeViewers )
{
    try { RCADM_d(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c* AHandle, bool includeViewers )
{
    try { RCADM_c(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
}

void ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z* AHandle, bool includeViewers )
{
    try { RCADM_z(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
}

// void DistMatrix<T,U,V>::Align( Int colAlign, Int rowAlign, bool constrain )
// ---------------------------------------------------------------------------
void ElDistMatrixAlign_s
( ElDistMatrix_s* AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { RCADM_s(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlign_d
( ElDistMatrix_d* AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { RCADM_d(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlign_c
( ElDistMatrix_c* AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { RCADM_c(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlign_z
( ElDistMatrix_z* AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { RCADM_z(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
}

// void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain )
// -----------------------------------------------------------------
void ElDistMatrixAlignCols_s
( ElDistMatrix_s* AHandle, ElInt colAlign, bool constrain )
{
    try { RCADM_s(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignCols_d
( ElDistMatrix_d* AHandle, ElInt colAlign, bool constrain )
{
    try { RCADM_d(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignCols_c
( ElDistMatrix_c* AHandle, ElInt colAlign, bool constrain )
{
    try { RCADM_c(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignCols_z
( ElDistMatrix_z* AHandle, ElInt colAlign, bool constrain )
{
    try { RCADM_z(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
}

// void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain )
// -----------------------------------------------------------------
void ElDistMatrixAlignRows_s
( ElDistMatrix_s* AHandle, ElInt rowAlign, bool constrain )
{
    try { RCADM_s(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignRows_d
( ElDistMatrix_d* AHandle, ElInt rowAlign, bool constrain )
{
    try { RCADM_d(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignRows_c
( ElDistMatrix_c* AHandle, ElInt rowAlign, bool constrain )
{
    try { RCADM_c(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
}

void ElDistMatrixAlignRows_z
( ElDistMatrix_z* AHandle, ElInt rowAlign, bool constrain )
{
    try { RCADM_z(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
}

// void DistMatrix<T,U,V>::FreeAlignments()
// ----------------------------------------
void ElDistMatrixFreeAlignments_s( ElDistMatrix_s* AHandle )
{ RCADM_s(AHandle)->FreeAlignments(); }

void ElDistMatrixFreeAlignments_d( ElDistMatrix_d* AHandle )
{ RCADM_d(AHandle)->FreeAlignments(); }

void ElDistMatrixFreeAlignments_c( ElDistMatrix_c* AHandle )
{ RCADM_c(AHandle)->FreeAlignments(); }

void ElDistMatrixFreeAlignments_z( ElDistMatrix_z* AHandle )
{ RCADM_z(AHandle)->FreeAlignments(); }

// void DistMatrix<T,U,V>::SetRoot( Int root )
// -------------------------------------------
void ElDistMatrixSetRoot_s( ElDistMatrix_s* AHandle, ElInt root )
{
    try { RCADM_s(AHandle)->SetRoot(root); }
    CATCH
}

void ElDistMatrixSetRoot_d( ElDistMatrix_d* AHandle, ElInt root )
{
    try { RCADM_d(AHandle)->SetRoot(root); }
    CATCH
}

void ElDistMatrixSetRoot_c( ElDistMatrix_c* AHandle, ElInt root )
{
    try { RCADM_c(AHandle)->SetRoot(root); }
    CATCH
}

void ElDistMatrixSetRoot_z( ElDistMatrix_z* AHandle, ElInt root )
{
    try { RCADM_z(AHandle)->SetRoot(root); }
    CATCH
}

// TODO: Align[Cols,Rows]With. Need a C version of DistData

// TODO: Align[Cols,Rows]AndResize

// void DistMatrix<T,U,V>::Attach
// ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
//   T* buffer, Int ldim, Int root )
// ----------------------------------------------------------------------
void ElDistMatrixAttach_s
( ElDistMatrix_s* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, float* buffer, ElInt ldim, ElInt root )
{
    try { RCADM_s(AHandle)->Attach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
}

void ElDistMatrixAttach_d
( ElDistMatrix_d* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, double* buffer, ElInt ldim, ElInt root )
{
    try { RCADM_d(AHandle)->Attach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
}

void ElDistMatrixAttach_c
( ElDistMatrix_c* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, complex_float* buffer, ElInt ldim, 
  ElInt root )
{
    try { RCADM_c(AHandle)->Attach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,RCB_c(buffer),
           ldim,root); }
    CATCH
}

void ElDistMatrixAttach_z
( ElDistMatrix_z* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, complex_double* buffer, ElInt ldim, 
  ElInt root )
{
    try { RCADM_z(AHandle)->Attach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,RCB_z(buffer),
           ldim,root); }
    CATCH
}

// void DistMatrix<T,U,V>::LockedAttach
// ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
//   const T* buffer, Int ldim, Int root )
// ----------------------------------------------------------------------
void ElDistMatrixLockedAttach_s
( ElDistMatrix_s* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const float* buffer, 
  ElInt ldim, ElInt root )
{
    try { RCADM_s(AHandle)->LockedAttach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
}

void ElDistMatrixLockedAttach_d
( ElDistMatrix_d* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const double* buffer, 
  ElInt ldim, ElInt root )
{
    try { RCADM_d(AHandle)->LockedAttach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
}

void ElDistMatrixLockedAttach_c
( ElDistMatrix_c* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const complex_float* buffer, 
  ElInt ldim, ElInt root )
{
    try { RCADM_c(AHandle)->LockedAttach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,
           RCB_c_const(buffer),ldim,root); }
    CATCH
}

void ElDistMatrixLockedAttach_z
( ElDistMatrix_z* AHandle, ElInt height, ElInt width, const ElGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const complex_double* buffer, 
  ElInt ldim, ElInt root )
{
    try { RCADM_z(AHandle)->LockedAttach
          (height,width,*RCG_const(gridHandle),colAlign,rowAlign,
           RCB_z_const(buffer),ldim,root); }
    CATCH
}

// Int DistMatrix<T,U,V>::Height() const
// -------------------------------------
ElInt ElDistMatrixHeight_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->Height(); }

ElInt ElDistMatrixHeight_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->Height(); }

ElInt ElDistMatrixHeight_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->Height(); }

ElInt ElDistMatrixHeight_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->Height(); }

// Int DistMatrix<T,U,V>::Width() const
// ------------------------------------
ElInt ElDistMatrixWidth_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->Width(); }

ElInt ElDistMatrixWidth_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->Width(); }

ElInt ElDistMatrixWidth_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->Width(); }

ElInt ElDistMatrixWidth_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->Width(); }

// Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const
// ---------------------------------------------------------
ElInt ElDistMatrixDiagonalLength_s
( const ElDistMatrix_s* AHandle, ElInt offset )
{ return RCADM_s_const(AHandle)->DiagonalLength(offset); }

ElInt ElDistMatrixDiagonalLength_d
( const ElDistMatrix_d* AHandle, ElInt offset )
{ return RCADM_d_const(AHandle)->DiagonalLength(offset); }

ElInt ElDistMatrixDiagonalLength_c
( const ElDistMatrix_c* AHandle, ElInt offset )
{ return RCADM_c_const(AHandle)->DiagonalLength(offset); }

ElInt ElDistMatrixDiagonalLength_z
( const ElDistMatrix_z* AHandle, ElInt offset )
{ return RCADM_z_const(AHandle)->DiagonalLength(offset); }

// bool DistMatrix<T,U,V>::Viewing() const
// ---------------------------------------
bool ElDistMatrixViewing_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->Viewing(); }

bool ElDistMatrixViewing_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->Viewing(); }

bool ElDistMatrixViewing_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->Viewing(); }

bool ElDistMatrixViewing_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->Viewing(); }

// bool DistMatrix<T,U,V>::Locked() const
// --------------------------------------
bool ElDistMatrixLocked_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->Locked(); }

bool ElDistMatrixLocked_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->Locked(); }

bool ElDistMatrixLocked_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->Locked(); }

bool ElDistMatrixLocked_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->Locked(); }

// Int DistMatrix<T,U,V>::LocalHeight() const
// ------------------------------------------
ElInt ElDistMatrixLocalHeight_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->LocalHeight(); }

ElInt ElDistMatrixLocalHeight_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->LocalHeight(); }

ElInt ElDistMatrixLocalHeight_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->LocalHeight(); }

ElInt ElDistMatrixLocalHeight_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->LocalHeight(); }

// Int DistMatrix<T,U,V>::LocalWidth() const
// -----------------------------------------
ElInt ElDistMatrixLocalWidth_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->LocalWidth(); }

ElInt ElDistMatrixLocalWidth_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->LocalWidth(); }

ElInt ElDistMatrixLocalWidth_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->LocalWidth(); }

ElInt ElDistMatrixLocalWidth_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->LocalWidth(); }

// Int DistMatrix<T,U,V>::LDim() const
// -----------------------------------------
ElInt ElDistMatrixLDim_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->LDim(); }

ElInt ElDistMatrixLDim_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->LDim(); }

ElInt ElDistMatrixLDim_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->LDim(); }

ElInt ElDistMatrixLDim_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->LDim(); }

// Matrix<T>& DistMatrix<T,U,V>::Matrix()
// --------------------------------------
ElMatrix_s* ElDistMatrixMatrix_s( ElDistMatrix_s* AHandle )
{
    ElMatrix_s* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<ElMatrix_s*>
                       (&RCADM_s(AHandle)->Matrix()); }
    CATCH
    return ALocHandle;
}

ElMatrix_d* ElDistMatrixMatrix_d( ElDistMatrix_d* AHandle )
{
    ElMatrix_d* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<ElMatrix_d*>
                       (&RCADM_d(AHandle)->Matrix()); }
    CATCH
    return ALocHandle;
}

ElMatrix_c* ElDistMatrixMatrix_c( ElDistMatrix_c* AHandle )
{
    ElMatrix_c* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<ElMatrix_c*>
                       (&RCADM_c(AHandle)->Matrix()); }
    CATCH
    return ALocHandle;
}

ElMatrix_z* ElDistMatrixMatrix_z( ElDistMatrix_z* AHandle )
{
    ElMatrix_z* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<ElMatrix_z*>
                       (&RCADM_z(AHandle)->Matrix()); }
    CATCH
    return ALocHandle;
}

// const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const
// --------------------------------------------------------
const ElMatrix_s* ElDistMatrixLockedMatrix_s( const ElDistMatrix_s* AHandle )
{
    const ElMatrix_s* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<const ElMatrix_s*>
                       (&RCADM_s_const(AHandle)->LockedMatrix()); }
    CATCH
    return ALocHandle;
}

const ElMatrix_d* ElDistMatrixLockedMatrix_d( const ElDistMatrix_d* AHandle )
{
    const ElMatrix_d* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<const ElMatrix_d*>
                       (&RCADM_d_const(AHandle)->LockedMatrix()); }
    CATCH
    return ALocHandle;
}

const ElMatrix_c* ElDistMatrixLockedMatrix_c( const ElDistMatrix_c* AHandle )
{
    const ElMatrix_c* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<const ElMatrix_c*>
                       (&RCADM_c_const(AHandle)->LockedMatrix()); }
    CATCH
    return ALocHandle;
}

const ElMatrix_z* ElDistMatrixLockedMatrix_z( const ElDistMatrix_z* AHandle )
{
    const ElMatrix_z* ALocHandle = 0;
    try { ALocHandle = reinterpret_cast<const ElMatrix_z*>
                       (&RCADM_z_const(AHandle)->LockedMatrix()); }
    CATCH
    return ALocHandle;
}

// size_t DistMatrix<T,U,V>::AllocatedMemory() const
// -------------------------------------------------
size_t ElDistMatrixAllocatedMemory_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->AllocatedMemory(); }

size_t ElDistMatrixAllocatedMemory_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->AllocatedMemory(); }

size_t ElDistMatrixAllocatedMemory_c( const ElDistMatrix_c* AHandle )
{ return RCADM_c_const(AHandle)->AllocatedMemory(); }

size_t ElDistMatrixAllocatedMemory_z( const ElDistMatrix_z* AHandle )
{ return RCADM_z_const(AHandle)->AllocatedMemory(); }

// T* DistMatrix<T,U,V>::Buffer()
// ------------------------------
float* ElDistMatrixBuffer_s( ElDistMatrix_s* AHandle )
{ return RCADM_s(AHandle)->Buffer(); }

double* ElDistMatrixBuffer_d( ElDistMatrix_d* AHandle )
{ return RCADM_d(AHandle)->Buffer(); }

complex_float* ElDistMatrixBuffer_c( ElDistMatrix_c* AHandle )
{ return (complex_float*)RCADM_c(AHandle)->Buffer(); }

complex_double* ElDistMatrixBuffer_z( ElDistMatrix_z* AHandle )
{ return (complex_double*)RCADM_z(AHandle)->Buffer(); }

// const T* DistMatrix<T,U,V>::LockedBuffer() const
// ------------------------------------------------
const float* 
ElDistMatrixLockedBuffer_s( const ElDistMatrix_s* AHandle )
{ return RCADM_s_const(AHandle)->LockedBuffer(); }

const double* 
ElDistMatrixLockedBuffer_d( const ElDistMatrix_d* AHandle )
{ return RCADM_d_const(AHandle)->LockedBuffer(); }

const complex_float* 
ElDistMatrixLockedBuffer_c( const ElDistMatrix_c* AHandle )
{ return (const complex_float*)RCADM_c_const(AHandle)->LockedBuffer(); }

const complex_double* 
ElDistMatrixLockedBuffer_z( const ElDistMatrix_z* AHandle )
{ return (const complex_double*)RCADM_z_const(AHandle)->LockedBuffer(); }

// const Grid& DistMatrix<T,U,V>::Grid() const
// -------------------------------------------
const ElGrid ElDistMatrixGrid_s( const ElDistMatrix_s* AHandle )
{ return (const ElGrid)reinterpret_cast<const struct ElGridDummy*>
         (&RCADM_s_const(AHandle)->Grid()); }

const ElGrid ElDistMatrixGrid_d( const ElDistMatrix_d* AHandle )
{ return (const ElGrid)reinterpret_cast<const struct ElGridDummy*>
         (&RCADM_d_const(AHandle)->Grid()); }

const ElGrid ElDistMatrixGrid_c( const ElDistMatrix_c* AHandle )
{ return (const ElGrid)reinterpret_cast<const struct ElGridDummy*>
         (&RCADM_c_const(AHandle)->Grid()); }

const ElGrid ElDistMatrixGrid_z( const ElDistMatrix_z* AHandle )
{ return (const ElGrid)reinterpret_cast<const struct ElGridDummy*>
         (&RCADM_z_const(AHandle)->Grid()); }

// void DistMatrix<T,U,V>::Get( Int i, Int j ) const
// -------------------------------------------------
float ElDistMatrixGet_s( const ElDistMatrix_s* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCADM_s_const(AHandle)->Get(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGet_d( const ElDistMatrix_d* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCADM_d_const(AHandle)->Get(i,j); }
    CATCH
    return alpha;
}

complex_float ElDistMatrixGet_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    complex_float alpha;
    try 
    { 
        Complex<float> alphaC = RCADM_c_const(AHandle)->Get(i,j);
        alpha.real = alphaC.real();
        alpha.imag = alphaC.imag();
    }
    CATCH
    return alpha;
}

complex_double ElDistMatrixGet_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    complex_double alpha;
    try 
    { 
        Complex<double> alphaC = RCADM_z_const(AHandle)->Get(i,j); 
        alpha.real = alphaC.real();
        alpha.imag = alphaC.imag();
    }
    CATCH
    return alpha;
}

// void DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
// ---------------------------------------------------------
float ElDistMatrixGetRealPart_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCADM_c_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGetRealPart_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCADM_c_const(AHandle)->GetRealPart(i,j); }
    CATCH
    return alpha;
}

// void DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
// ---------------------------------------------------------
float ElDistMatrixGetImagPart_c
( const ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    float alpha = -1;
    try { alpha = RCADM_c_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return alpha;
}

double ElDistMatrixGetImagPart_z
( const ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    double alpha = -1;
    try { alpha = RCADM_c_const(AHandle)->GetImagPart(i,j); }
    CATCH
    return alpha;
}

// void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
// ----------------------------------------------------
void ElDistMatrixSet_s( ElDistMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_s(AHandle)->Set(i,j,alpha); }
    CATCH
}

void ElDistMatrixSet_d
( ElDistMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_d(AHandle)->Set(i,j,alpha); }
    CATCH
}

void ElDistMatrixSet_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { RCADM_c(AHandle)->Set(i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
}

void ElDistMatrixSet_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { RCADM_z(AHandle)->Set(i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
}

// void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetRealPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_c(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixSetRealPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_z(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
void ElDistMatrixSetImagPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_c(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixSetImagPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_z(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
// -------------------------------------------------------
void ElDistMatrixUpdate_s
( ElDistMatrix_s* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_s(AHandle)->Update(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdate_d
( ElDistMatrix_d* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_d(AHandle)->Update(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdate_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { RCADM_c(AHandle)->Update
          (i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
}

void ElDistMatrixUpdate_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { RCADM_z(AHandle)->Update
          (i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
}

// void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_c(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_z(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
void ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c* AHandle, ElInt i, ElInt j, float alpha )
{
    try { RCADM_c(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
}

void ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z* AHandle, ElInt i, ElInt j, double alpha )
{
    try { RCADM_z(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
}

// void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
// ------------------------------------------------
void ElDistMatrixMakeReal_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    try { RCADM_c(AHandle)->MakeReal(i,j); }
    CATCH
}

void ElDistMatrixMakeReal_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    try { RCADM_z(AHandle)->MakeReal(i,j); }
    CATCH
}

// void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
// -------------------------------------------------
void ElDistMatrixConjugate_c( ElDistMatrix_c* AHandle, ElInt i, ElInt j )
{
    try { RCADM_c(AHandle)->Conjugate(i,j); }
    CATCH
}

void ElDistMatrixConjugate_z( ElDistMatrix_z* AHandle, ElInt i, ElInt j )
{
    try { RCADM_z(AHandle)->Conjugate(i,j); }
    CATCH
}

// DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
// ----------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetDiagonal_s
( const ElDistMatrix_s* AHandle, ElInt offset )
{
    ElDistMatrix_s* dHandle = 0;
    try 
    {
        auto AAbs = RCADM_s_const(AHandle);
        const Grid& grid = AAbs->Grid();

        const DistData distData = AAbs->DistData();
        const Dist U = distData.colDist;
        const Dist V = distData.rowDist;

        if( U == CIRC && V == CIRC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,CIRC,CIRC>*>(AAbs);
            auto* d = new DistMatrix<float,CIRC,CIRC>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == MC && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,MC,MR>*>(AAbs);
            auto* d = new DistMatrix<float,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == MC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,MC,STAR>*>(AAbs);
            auto* d = new DistMatrix<float,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == MD && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,MD,STAR>*>(AAbs);
            auto* d = new DistMatrix<float,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == MC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,MC>*>(AAbs);
            auto* d = new DistMatrix<float,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == MD )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,MD>*>(AAbs);
            auto* d = new DistMatrix<float,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,MR>*>(AAbs);
            auto* d = new DistMatrix<float,MR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,STAR>*>(AAbs);
            auto* d = new DistMatrix<float,STAR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == VC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,VC>*>(AAbs);
            auto* d = new DistMatrix<float,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == STAR && V == VR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,STAR,VR>*>(AAbs);
            auto* d = new DistMatrix<float,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == VC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,VC,STAR>*>(AAbs);
            auto* d = new DistMatrix<float,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else if( U == VR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<float,VR,STAR>*>(AAbs);
            auto* d = new DistMatrix<float,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_s*>(d);
        }
        else
            RuntimeError("Invalid distribution pair");
    }
    CATCH
    return dHandle;
}

ElDistMatrix_d* ElDistMatrixGetDiagonal_d
( const ElDistMatrix_d* AHandle, ElInt offset )
{
    ElDistMatrix_d* dHandle = 0;
    try 
    {
        auto AAbs = RCADM_d_const(AHandle);
        const Grid& grid = AAbs->Grid();

        const DistData distData = AAbs->DistData();
        const Dist U = distData.colDist;
        const Dist V = distData.rowDist;

        if( U == CIRC && V == CIRC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,CIRC,CIRC>*>(AAbs);
            auto* d = new DistMatrix<double,CIRC,CIRC>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == MC && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,MC,MR>*>(AAbs);
            auto* d = new DistMatrix<double,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == MC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,MC,STAR>*>(AAbs);
            auto* d = new DistMatrix<double,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == MD && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,MD,STAR>*>(AAbs);
            auto* d = new DistMatrix<double,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == MC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,MC>*>(AAbs);
            auto* d = new DistMatrix<double,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == MD )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,MD>*>(AAbs);
            auto* d = new DistMatrix<double,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,MR>*>(AAbs);
            auto* d = new DistMatrix<double,MR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,STAR>*>(AAbs);
            auto* d = new DistMatrix<double,STAR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == VC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,VC>*>(AAbs);
            auto* d = new DistMatrix<double,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == STAR && V == VR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,STAR,VR>*>(AAbs);
            auto* d = new DistMatrix<double,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == VC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,VC,STAR>*>(AAbs);
            auto* d = new DistMatrix<double,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else if( U == VR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<double,VR,STAR>*>(AAbs);
            auto* d = new DistMatrix<double,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_d*>(d);
        }
        else
            RuntimeError("Invalid distribution pair");
    }
    CATCH
    return dHandle;
}

ElDistMatrix_c* ElDistMatrixGetDiagonal_c
( const ElDistMatrix_c* AHandle, ElInt offset )
{
    ElDistMatrix_c* dHandle = 0;
    try 
    {
        auto AAbs = RCADM_c_const(AHandle);
        const Grid& grid = AAbs->Grid();

        const DistData distData = AAbs->DistData();
        const Dist U = distData.colDist;
        const Dist V = distData.rowDist;

        if( U == CIRC && V == CIRC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,CIRC,CIRC>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,CIRC,CIRC>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == MC && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,MC,MR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == MC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,MC,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == MD && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,MD,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == MC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,MC>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == MD )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,MD>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,MR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,MR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,STAR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == VC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,VC>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == STAR && V == VR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,STAR,VR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == VC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,VC,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else if( U == VR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<float>,VR,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<float>,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_c*>(d);
        }
        else
            RuntimeError("Invalid distribution pair");
    }
    CATCH
    return dHandle;
}

ElDistMatrix_z* ElDistMatrixGetDiagonal_z
( const ElDistMatrix_z* AHandle, ElInt offset )
{
    ElDistMatrix_z* dHandle = 0;
    try 
    {
        auto AAbs = RCADM_z_const(AHandle);
        const Grid& grid = AAbs->Grid();

        const DistData distData = AAbs->DistData();
        const Dist U = distData.colDist;
        const Dist V = distData.rowDist;

        if( U == CIRC && V == CIRC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,CIRC,CIRC>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,CIRC,CIRC>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == MC && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,MC,MR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == MC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,MC,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == MD && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,MD,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == MC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,MC>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == MD )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,MD>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MD,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == MR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,MR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,MR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,STAR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == VC )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,VC>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == STAR && V == VR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,STAR,VR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == VC && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,VC,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,VC,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else if( U == VR && V == STAR )
        {
            auto A = 
              dynamic_cast<const DistMatrix<Complex<double>,VR,STAR>*>(AAbs);
            auto* d = new DistMatrix<Complex<double>,VR,STAR>(grid);
            A->GetDiagonal( *d, offset );
            dHandle = reinterpret_cast<ElDistMatrix_z*>(d);
        }
        else
            RuntimeError("Invalid distribution pair");
    }
    CATCH
    return dHandle;
}

// TODO: More diagonal manipulation
// ================================

// DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElDistMatrix_s* ElDistMatrixGetSubmatrix_s
( const ElDistMatrix_s* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_s* ASubHandle = 0;
    try
    {
        const ElGrid gridHandle = ElDistMatrixGrid_s( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_s
                     ( EL_STAR, EL_STAR, gridHandle );

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCADM_s_const(AHandle)->GetSubmatrix
        ( rowIndVec, colIndVec, 
          *reinterpret_cast<DistMatrix<float,STAR,STAR>*>(ASubHandle) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_d* ElDistMatrixGetSubmatrix_d
( const ElDistMatrix_d* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_d* ASubHandle = 0;
    try
    {
        const ElGrid gridHandle = ElDistMatrixGrid_d( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_d
                     ( EL_STAR, EL_STAR, gridHandle );

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCADM_d_const(AHandle)->GetSubmatrix
        ( rowIndVec, colIndVec,
          *reinterpret_cast<DistMatrix<double,STAR,STAR>*>(ASubHandle) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_c* ElDistMatrixGetSubmatrix_c
( const ElDistMatrix_c* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_c* ASubHandle = 0;
    try
    {
        const ElGrid gridHandle = ElDistMatrixGrid_c( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_c
                     ( EL_STAR, EL_STAR, gridHandle );

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCADM_c_const(AHandle)->GetSubmatrix
        ( rowIndVec, colIndVec,
          *reinterpret_cast<DistMatrix<Complex<float>,STAR,STAR>*>
           (ASubHandle) );
    }
    CATCH
    return ASubHandle;
}

ElDistMatrix_z* ElDistMatrixGetSubmatrix_z
( const ElDistMatrix_z* AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds )
{
    ElDistMatrix_z* ASubHandle = 0;
    try
    {
        const ElGrid gridHandle = ElDistMatrixGrid_z( AHandle );
        ASubHandle = ElDistMatrixCreateSpecific_z
                     ( EL_STAR, EL_STAR, gridHandle );

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        RCADM_z_const(AHandle)->GetSubmatrix
        ( rowIndVec, colIndVec,
          *reinterpret_cast<DistMatrix<Complex<double>,STAR,STAR>*>
           (ASubHandle) );
    }
    CATCH
    return ASubHandle;
}

/// TODO: Many more member functions

} // extern "C"
