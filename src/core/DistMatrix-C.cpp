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
#include "El-C.h"
using namespace El;

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

extern "C" {

// DistMatrix<T,MC,MR>::DistMatrix( const Grid& g )
// ------------------------------------------------
ElError ElDistMatrixCreate_s( ElConstGrid gridHandle, ElDistMatrix_s* AHandle )
{
    try 
    {
        auto grid = Reinterpret(gridHandle);
        *AHandle = Reinterpret(new DistMatrix<float>(*grid));
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCreate_d( ElConstGrid gridHandle, ElDistMatrix_d* AHandle )
{
    try 
    {
        auto grid = Reinterpret(gridHandle);
        *AHandle = Reinterpret(new DistMatrix<double>(*grid));
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCreate_c( ElConstGrid gridHandle, ElDistMatrix_c* AHandle )
{
    try 
    {
        auto grid = Reinterpret(gridHandle);
        *AHandle = Reinterpret(new DistMatrix<Complex<float>>(*grid));
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCreate_z( ElConstGrid gridHandle, ElDistMatrix_z* AHandle )
{
    try 
    {
        auto grid = Reinterpret(gridHandle);
        *AHandle = Reinterpret(new DistMatrix<Complex<double>>(*grid));
    }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"

// DistMatrix<T,U,V>::DistMatrix( const Grid& g )
// ----------------------------------------------

template<typename T>
ElError ElDistMatrixCreateSpecific
( ElDist U_C, ElDist V_C, ElConstGrid gridHandle, AbstractDistMatrix<T>** A )
{
    try 
    {
        Dist U = Reinterpret(U_C);
        Dist V = Reinterpret(V_C);
        auto grid = Reinterpret(gridHandle);

        if( U == CIRC && V == CIRC )
            *A = new DistMatrix<T,CIRC,CIRC>(*grid); 
        else if( U == MC && V == MR )
            *A = new DistMatrix<T,MC,MR>(*grid);
        else if( U == MC && V == STAR )
            *A = new DistMatrix<T,MC,STAR>(*grid);
        else if( U == MD && V == STAR )
            *A = new DistMatrix<T,MD,STAR>(*grid);
        else if( U == MR && V == MC )
            *A = new DistMatrix<T,MR,MC>(*grid);
        else if( U == MR && V == STAR )
            *A = new DistMatrix<T,MR,STAR>(*grid);
        else if( U == STAR && V == MC )
            *A = new DistMatrix<T,STAR,MC>(*grid);
        else if( U == STAR && V == MD )
            *A = new DistMatrix<T,STAR,MD>(*grid);
        else if( U == STAR && V == MR )
            *A = new DistMatrix<T,STAR,MR>(*grid);
        else if( U == STAR && V == STAR )
            *A = new DistMatrix<T,STAR,STAR>(*grid);
        else if( U == STAR && V == VC )
            *A = new DistMatrix<T,STAR,VC>(*grid);
        else if( U == STAR && V == VR )
            *A = new DistMatrix<T,STAR,VR>(*grid);
        else if( U == VC && V == STAR )
            *A = new DistMatrix<T,VC,STAR>(*grid);
        else if( U == VR && V == STAR )
            *A = new DistMatrix<T,VR,STAR>(*grid);
    }
    CATCH
    return EL_SUCCESS;
}

extern "C" {

ElError ElDistMatrixCreateSpecific_s
( ElDist U, ElDist V, ElConstGrid gridHandle, ElDistMatrix_s* AHandle )
{
    AbstractDistMatrix<float>* ADM;
    ElError error = ElDistMatrixCreateSpecific( U, V, gridHandle, &ADM );
    *AHandle = Reinterpret(ADM);
    return error;
}

ElError ElDistMatrixCreateSpecific_d
( ElDist U, ElDist V, ElConstGrid gridHandle, ElDistMatrix_d* AHandle )
{
    AbstractDistMatrix<double>* ADM;
    ElError error = ElDistMatrixCreateSpecific( U, V, gridHandle, &ADM );
    *AHandle = Reinterpret(ADM);
    return error;
}

ElError ElDistMatrixCreateSpecific_c
( ElDist U, ElDist V, ElConstGrid gridHandle, ElDistMatrix_c* AHandle )
{
    AbstractDistMatrix<Complex<float>>* ADM;
    ElError error = ElDistMatrixCreateSpecific( U, V, gridHandle, &ADM );
    *AHandle = Reinterpret(ADM);
    return error;
}

ElError ElDistMatrixCreateSpecific_z
( ElDist U, ElDist V, ElConstGrid gridHandle, ElDistMatrix_z* AHandle )
{
    AbstractDistMatrix<Complex<double>>* ADM;
    ElError error = ElDistMatrixCreateSpecific( U, V, gridHandle, &ADM );
    *AHandle = Reinterpret(ADM);
    return error;
}

// DistMatrix<T,U,V>::~DistMatrix()
// --------------------------------
ElError ElDistMatrixDestroy_s( ElConstDistMatrix_s AHandle )
{ 
    try { delete Reinterpret(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDestroy_d( ElConstDistMatrix_d AHandle )
{ 
    try { delete Reinterpret(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDestroy_c( ElConstDistMatrix_c AHandle )
{ 
    try { delete Reinterpret(AHandle); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDestroy_z( ElConstDistMatrix_z AHandle )
{ 
    try { delete Reinterpret(AHandle); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Empty()
// -------------------------------
ElError ElDistMatrixEmpty_s( ElDistMatrix_s AHandle )
{
    try { Reinterpret(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmpty_d( ElDistMatrix_d AHandle )
{
    try { Reinterpret(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmpty_c( ElDistMatrix_c AHandle )
{
    try { Reinterpret(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmpty_z( ElDistMatrix_z AHandle )
{
    try { Reinterpret(AHandle)->Empty(); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::EmptyData()
// -----------------------------------
ElError ElDistMatrixEmptyData_s( ElDistMatrix_s AHandle )
{
    try { Reinterpret(AHandle)->EmptyData(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmptyData_d( ElDistMatrix_d AHandle )
{
    try { Reinterpret(AHandle)->EmptyData(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmptyData_c( ElDistMatrix_c AHandle )
{
    try { Reinterpret(AHandle)->EmptyData(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixEmptyData_z( ElDistMatrix_z AHandle )
{
    try { Reinterpret(AHandle)->EmptyData(); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetGrid( const Grid& g )
// ------------------------------------------------
ElError ElDistMatrixSetGrid_s( ElDistMatrix_s AHandle, ElConstGrid gridHandle )
{
    try { Reinterpret(AHandle)->SetGrid(*Reinterpret(gridHandle)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetGrid_d( ElDistMatrix_d AHandle, ElConstGrid gridHandle )
{
    try { Reinterpret(AHandle)->SetGrid(*Reinterpret(gridHandle)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetGrid_c( ElDistMatrix_c AHandle, ElConstGrid gridHandle )
{
    try { Reinterpret(AHandle)->SetGrid(*Reinterpret(gridHandle)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetGrid_z( ElDistMatrix_z AHandle, ElConstGrid gridHandle )
{
    try { Reinterpret(AHandle)->SetGrid(*Reinterpret(gridHandle)); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Resize( Int height, Int width )
// -------------------------------------------------------
ElError ElDistMatrixResize_s
( ElDistMatrix_s AHandle, ElInt height, ElInt width )
{
    try { Reinterpret(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResize_d
( ElDistMatrix_d AHandle, ElInt height, ElInt width )
{
    try { Reinterpret(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResize_c
( ElDistMatrix_c AHandle, ElInt height, ElInt width )
{
    try { Reinterpret(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResize_z
( ElDistMatrix_z AHandle, ElInt height, ElInt width )
{
    try { Reinterpret(AHandle)->Resize(height,width); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Resize( Int height, Int width, Int ldim )
// -----------------------------------------------------------------
ElError ElDistMatrixResizeWithLDim_s
( ElDistMatrix_s AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { Reinterpret(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResizeWithLDim_d
( ElDistMatrix_d AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { Reinterpret(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResizeWithLDim_c
( ElDistMatrix_c AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { Reinterpret(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixResizeWithLDim_z
( ElDistMatrix_z AHandle, ElInt height, ElInt width, ElInt ldim )
{
    try { Reinterpret(AHandle)->Resize(height,width,ldim); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeConsistent()
// ----------------------------------------
ElError ElDistMatrixMakeConsistent_s
( ElDistMatrix_s AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeConsistent_d
( ElDistMatrix_d AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeConsistent_c
( ElDistMatrix_c AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeConsistent_z
( ElDistMatrix_z AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeSizeConsistent()
// --------------------------------------------
ElError ElDistMatrixMakeSizeConsistent_s
( ElDistMatrix_s AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeSizeConsistent_d
( ElDistMatrix_d AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeSizeConsistent_c
( ElDistMatrix_c AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeSizeConsistent_z
( ElDistMatrix_z AHandle, bool includeViewers )
{
    try { Reinterpret(AHandle)->MakeSizeConsistent(includeViewers); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Align( Int colAlign, Int rowAlign, bool constrain )
// ---------------------------------------------------------------------------
ElError ElDistMatrixAlign_s
( ElDistMatrix_s AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlign_d
( ElDistMatrix_d AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlign_c
( ElDistMatrix_c AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlign_z
( ElDistMatrix_z AHandle, ElInt colAlign, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->Align(colAlign,rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignCols( Int colAlign, bool constrain )
// -----------------------------------------------------------------
ElError ElDistMatrixAlignCols_s
( ElDistMatrix_s AHandle, ElInt colAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignCols_d
( ElDistMatrix_d AHandle, ElInt colAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignCols_c
( ElDistMatrix_c AHandle, ElInt colAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignCols_z
( ElDistMatrix_z AHandle, ElInt colAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignCols(colAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignRows( Int rowAlign, bool constrain )
// -----------------------------------------------------------------
ElError ElDistMatrixAlignRows_s
( ElDistMatrix_s AHandle, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRows_d
( ElDistMatrix_d AHandle, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRows_c
( ElDistMatrix_c AHandle, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRows_z
( ElDistMatrix_z AHandle, ElInt rowAlign, bool constrain )
{
    try { Reinterpret(AHandle)->AlignRows(rowAlign,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::FreeAlignments()
// ----------------------------------------
ElError ElDistMatrixFreeAlignments_s( ElDistMatrix_s AHandle )
{ 
    try { Reinterpret(AHandle)->FreeAlignments(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixFreeAlignments_d( ElDistMatrix_d AHandle )
{ 
    try { Reinterpret(AHandle)->FreeAlignments(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixFreeAlignments_c( ElDistMatrix_c AHandle )
{ 
    try { Reinterpret(AHandle)->FreeAlignments(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixFreeAlignments_z( ElDistMatrix_z AHandle )
{ 
    try { Reinterpret(AHandle)->FreeAlignments(); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetRoot( Int root )
// -------------------------------------------
ElError ElDistMatrixSetRoot_s( ElDistMatrix_s AHandle, ElInt root )
{
    try { Reinterpret(AHandle)->SetRoot(root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRoot_d( ElDistMatrix_d AHandle, ElInt root )
{
    try { Reinterpret(AHandle)->SetRoot(root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRoot_c( ElDistMatrix_c AHandle, ElInt root )
{
    try { Reinterpret(AHandle)->SetRoot(root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRoot_z( ElDistMatrix_z AHandle, ElInt root )
{
    try { Reinterpret(AHandle)->SetRoot(root); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignWith( const DistData& data, bool constrain )
// -------------------------------------------------------------------------
ElError ElDistMatrixAlignWith_s
( ElDistMatrix_s AHandle, ElDistData distData, bool constrain )
{
    try { Reinterpret(AHandle)->AlignWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignWith_d
( ElDistMatrix_d AHandle, ElDistData distData, bool constrain )
{
    try { Reinterpret(AHandle)->AlignWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignWith_c
( ElDistMatrix_c AHandle, ElDistData distData, bool constrain )
{
    try { Reinterpret(AHandle)->AlignWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignWith_z
( ElDistMatrix_z AHandle, ElDistData distData, bool constrain )
{
    try { Reinterpret(AHandle)->AlignWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignColsWith( const DistData& data, bool constrain )
// -----------------------------------------------------------------------------
ElError ElDistMatrixAlignColsWith_s
( ElDistMatrix_s AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsWith_d
( ElDistMatrix_d AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsWith_c
( ElDistMatrix_c AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsWith_z
( ElDistMatrix_z AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignRowsWith( const DistData& data, bool constrain )
// -----------------------------------------------------------------------------
ElError ElDistMatrixAlignRowsWith_s
( ElDistMatrix_s AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsWith_d
( ElDistMatrix_d AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsWith_c
( ElDistMatrix_c AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsWith_z
( ElDistMatrix_z AHandle, ElDistData distData, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsWith( Reinterpret(distData), constrain ); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignAndResize
// ( Int colAlign, Int rowAlign, Int height, Int width, 
//   bool force, bool constrain )
// ----------------------------------------------------
ElError ElDistMatrixAlignAndResize_s
( ElDistMatrix_s AHandle, 
  ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignAndResize
      (colAlign,rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignAndResize_d
( ElDistMatrix_d AHandle, 
  ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignAndResize
      (colAlign,rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignAndResize_c
( ElDistMatrix_c AHandle, 
  ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignAndResize
      (colAlign,rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignAndResize_z
( ElDistMatrix_z AHandle, 
  ElInt colAlign, ElInt rowAlign, ElInt height, ElInt width, 
  bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignAndResize
      (colAlign,rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignColsAndResize
// ( Int colAlign, Int height, Int width, bool force, bool constrain )
// -------------------------------------------------------------------
ElError ElDistMatrixAlignColsAndResize_s
( ElDistMatrix_s AHandle, 
  ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsAndResize
      (colAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsAndResize_d
( ElDistMatrix_d AHandle, 
  ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsAndResize
      (colAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsAndResize_c
( ElDistMatrix_c AHandle, 
  ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsAndResize
      (colAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignColsAndResize_z
( ElDistMatrix_z AHandle, 
  ElInt colAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignColsAndResize
      (colAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::AlignRowsAndResize
// ( Int rowAlign, Int height, Int width, bool force, bool constrain )
// -------------------------------------------------------------------
ElError ElDistMatrixAlignRowsAndResize_s
( ElDistMatrix_s AHandle, 
  ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsAndResize
      (rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsAndResize_d
( ElDistMatrix_d AHandle, 
  ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsAndResize
      (rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsAndResize_c
( ElDistMatrix_c AHandle, 
  ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsAndResize
      (rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAlignRowsAndResize_z
( ElDistMatrix_z AHandle, 
  ElInt rowAlign, ElInt height, ElInt width, bool force, bool constrain )
{
    try 
    { Reinterpret(AHandle)->AlignRowsAndResize
      (rowAlign,height,width,force,constrain); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Attach
// ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
//   T* buffer, Int ldim, Int root )
// ----------------------------------------------------------------------
ElError ElDistMatrixAttach_s
( ElDistMatrix_s AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, float* buffer, ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->Attach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAttach_d
( ElDistMatrix_d AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, double* buffer, ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->Attach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAttach_c
( ElDistMatrix_c AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, complex_float* buffer, ElInt ldim, 
  ElInt root )
{
    try { Reinterpret(AHandle)->Attach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,
           Reinterpret(buffer),ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAttach_z
( ElDistMatrix_z AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, complex_double* buffer, ElInt ldim, 
  ElInt root )
{
    try { Reinterpret(AHandle)->Attach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,
           Reinterpret(buffer),ldim,root); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::LockedAttach
// ( Int height, Int width, const Grid& grid, Int colAlign, Int rowAlign, 
//   const T* buffer, Int ldim, Int root )
// ----------------------------------------------------------------------
ElError ElDistMatrixLockedAttach_s
( ElDistMatrix_s AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const float* buffer, 
  ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->LockedAttach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedAttach_d
( ElDistMatrix_d AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const double* buffer, 
  ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->LockedAttach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,buffer,
           ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedAttach_c
( ElDistMatrix_c AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const complex_float* buffer, 
  ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->LockedAttach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,
           Reinterpret(buffer),ldim,root); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedAttach_z
( ElDistMatrix_z AHandle, ElInt height, ElInt width, ElConstGrid gridHandle,
  ElInt colAlign, ElInt rowAlign, const complex_double* buffer, 
  ElInt ldim, ElInt root )
{
    try { Reinterpret(AHandle)->LockedAttach
          (height,width,*Reinterpret(gridHandle),colAlign,rowAlign,
           Reinterpret(buffer),ldim,root); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::Height() const
// -------------------------------------
ElError ElDistMatrixHeight_s( ElConstDistMatrix_s AHandle, ElInt* height )
{ 
    try { *height = Reinterpret(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixHeight_d( ElConstDistMatrix_d AHandle, ElInt* height )
{ 
    try { *height = Reinterpret(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixHeight_c( ElConstDistMatrix_c AHandle, ElInt* height )
{ 
    try { *height = Reinterpret(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixHeight_z( ElConstDistMatrix_z AHandle, ElInt* height )
{ 
    try { *height = Reinterpret(AHandle)->Height(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::Width() const
// ------------------------------------
ElError ElDistMatrixWidth_s( ElConstDistMatrix_s AHandle, ElInt* width )
{ 
    try { *width = Reinterpret(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixWidth_d( ElConstDistMatrix_d AHandle, ElInt* width )
{ 
    try { *width = Reinterpret(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixWidth_c( ElConstDistMatrix_c AHandle, ElInt* width )
{ 
    try { *width = Reinterpret(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixWidth_z( ElConstDistMatrix_z AHandle, ElInt* width )
{ 
    try { *width = Reinterpret(AHandle)->Width(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::DiagonalLength( Int offset ) const
// ---------------------------------------------------------
ElError ElDistMatrixDiagonalLength_s
( ElConstDistMatrix_s AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = Reinterpret(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalLength_d
( ElConstDistMatrix_d AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = Reinterpret(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalLength_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = Reinterpret(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalLength_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElInt* length )
{ 
    try { *length = Reinterpret(AHandle)->DiagonalLength(offset); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::Viewing() const
// ---------------------------------------
ElError ElDistMatrixViewing_s( ElConstDistMatrix_s AHandle, bool* viewing )
{ 
    try { *viewing = Reinterpret(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixViewing_d( ElConstDistMatrix_d AHandle, bool* viewing )
{ 
    try { *viewing = Reinterpret(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixViewing_c( ElConstDistMatrix_c AHandle, bool* viewing )
{ 
    try { *viewing = Reinterpret(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixViewing_z( ElConstDistMatrix_z AHandle, bool* viewing )
{ 
    try { *viewing = Reinterpret(AHandle)->Viewing(); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::Locked() const
// --------------------------------------
ElError ElDistMatrixLocked_s( ElConstDistMatrix_s AHandle, bool* locked )
{ 
    try { *locked = Reinterpret(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocked_d( ElConstDistMatrix_d AHandle, bool* locked )
{ 
    try { *locked = Reinterpret(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocked_c( ElConstDistMatrix_c AHandle, bool* locked )
{ 
    try { *locked = Reinterpret(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocked_z( ElConstDistMatrix_z AHandle, bool* locked )
{ 
    try { *locked = Reinterpret(AHandle)->Locked(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalHeight() const
// ------------------------------------------
ElError ElDistMatrixLocalHeight_s
( ElConstDistMatrix_s AHandle, ElInt* localHeight )
{ 
    try { *localHeight = Reinterpret(AHandle)->LocalHeight(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalHeight_d
( ElConstDistMatrix_d AHandle, ElInt* localHeight )
{ 
    try { *localHeight = Reinterpret(AHandle)->LocalHeight(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalHeight_c
( ElConstDistMatrix_c AHandle, ElInt* localHeight )
{ 
    try { *localHeight = Reinterpret(AHandle)->LocalHeight(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalHeight_z
( ElConstDistMatrix_z AHandle, ElInt* localHeight )
{ 
    try { *localHeight = Reinterpret(AHandle)->LocalHeight(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalWidth() const
// -----------------------------------------
ElError ElDistMatrixLocalWidth_s
( ElConstDistMatrix_s AHandle, ElInt* localWidth )
{ 
    try { *localWidth = Reinterpret(AHandle)->LocalWidth(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalWidth_d
( ElConstDistMatrix_d AHandle, ElInt* localWidth )
{ 
    try { *localWidth = Reinterpret(AHandle)->LocalWidth(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalWidth_c
( ElConstDistMatrix_c AHandle, ElInt* localWidth )
{ 
    try { *localWidth = Reinterpret(AHandle)->LocalWidth(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalWidth_z
( ElConstDistMatrix_z AHandle, ElInt* localWidth )
{ 
    try { *localWidth = Reinterpret(AHandle)->LocalWidth(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LDim() const
// -----------------------------------------
ElError ElDistMatrixLDim_s( ElConstDistMatrix_s AHandle, ElInt* ldim )
{ 
    try { *ldim = Reinterpret(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS; 
}

ElError ElDistMatrixLDim_d( ElConstDistMatrix_d AHandle, ElInt* ldim )
{ 
    try { *ldim = Reinterpret(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLDim_c( ElConstDistMatrix_c AHandle, ElInt* ldim )
{ 
    try { *ldim = Reinterpret(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLDim_z( ElConstDistMatrix_z AHandle, ElInt* ldim )
{ 
    try { *ldim = Reinterpret(AHandle)->LDim(); }
    CATCH
    return EL_SUCCESS;
}

// Matrix<T>& DistMatrix<T,U,V>::Matrix()
// --------------------------------------
ElError ElDistMatrixMatrix_s( ElDistMatrix_s AHandle, ElMatrix_s* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->Matrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMatrix_d( ElDistMatrix_d AHandle, ElMatrix_d* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->Matrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMatrix_c( ElDistMatrix_c AHandle, ElMatrix_c* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->Matrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMatrix_z( ElDistMatrix_z AHandle, ElMatrix_z* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->Matrix()); }
    CATCH
    return EL_SUCCESS;
}

// const Matrix<T>& DistMatrix<T,U,V>::LockedMatrix() const
// --------------------------------------------------------
ElError ElDistMatrixLockedMatrix_s
( ElConstDistMatrix_s AHandle, ElConstMatrix_s* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->LockedMatrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedMatrix_d
( ElConstDistMatrix_d AHandle, ElConstMatrix_d* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->LockedMatrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedMatrix_c
( ElConstDistMatrix_c AHandle, ElConstMatrix_c* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->LockedMatrix()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedMatrix_z
( ElConstDistMatrix_z AHandle, ElConstMatrix_z* ALocHandle )
{
    try { *ALocHandle = Reinterpret(&Reinterpret(AHandle)->LockedMatrix()); }
    CATCH
    return EL_SUCCESS;
}

// size_t DistMatrix<T,U,V>::AllocatedMemory() const
// -------------------------------------------------
ElError ElDistMatrixAllocatedMemory_s
( ElConstDistMatrix_s AHandle, size_t* mem )
{ 
    try { *mem = Reinterpret(AHandle)->AllocatedMemory(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAllocatedMemory_d
( ElConstDistMatrix_d AHandle, size_t* mem )
{ 
    try { *mem = Reinterpret(AHandle)->AllocatedMemory(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAllocatedMemory_c
( ElConstDistMatrix_c AHandle, size_t* mem )
{ 
    try { *mem = Reinterpret(AHandle)->AllocatedMemory(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixAllocatedMemory_z
( ElConstDistMatrix_z AHandle, size_t* mem )
{ 
    try { *mem = Reinterpret(AHandle)->AllocatedMemory(); }
    CATCH
    return EL_SUCCESS;
}

// T* DistMatrix<T,U,V>::Buffer()
// ------------------------------
ElError ElDistMatrixBuffer_s( ElDistMatrix_s AHandle, float** buffer )
{ 
    try { *buffer = Reinterpret(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixBuffer_d( ElDistMatrix_d AHandle, double** buffer )
{ 
    try { *buffer = Reinterpret(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixBuffer_c( ElDistMatrix_c AHandle, complex_float** buffer )
{ 
    try { *buffer = (complex_float*)Reinterpret(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixBuffer_z( ElDistMatrix_z AHandle, complex_double** buffer )
{ 
    try { *buffer = (complex_double*)Reinterpret(AHandle)->Buffer(); }
    CATCH
    return EL_SUCCESS;
}

// const T* DistMatrix<T,U,V>::LockedBuffer() const
// ------------------------------------------------
ElError ElDistMatrixLockedBuffer_s
( ElConstDistMatrix_s AHandle, const float** buffer )
{ 
    try { *buffer = Reinterpret(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedBuffer_d
( ElConstDistMatrix_d AHandle, const double** buffer )
{ 
    try { *buffer = Reinterpret(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedBuffer_c
( ElConstDistMatrix_c AHandle, const complex_float** buffer )
{ 
    try 
    { *buffer = (const complex_float*)Reinterpret(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLockedBuffer_z
( ElConstDistMatrix_z AHandle, const complex_double** buffer )
{ 
    try 
    { *buffer = (const complex_double*)Reinterpret(AHandle)->LockedBuffer(); }
    CATCH
    return EL_SUCCESS;
}

// const Grid& DistMatrix<T,U,V>::Grid() const
// -------------------------------------------
ElError ElDistMatrixGrid_s
( ElConstDistMatrix_s AHandle, ElConstGrid* gridHandle )
{ 
    try { *gridHandle = Reinterpret(&Reinterpret(AHandle)->Grid()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGrid_d
( ElConstDistMatrix_d AHandle, ElConstGrid* gridHandle )
{ 
    try { *gridHandle = Reinterpret(&Reinterpret(AHandle)->Grid()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGrid_c
( ElConstDistMatrix_c AHandle, ElConstGrid* gridHandle )
{ 
    try { *gridHandle = Reinterpret(&Reinterpret(AHandle)->Grid()); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGrid_z
( ElConstDistMatrix_z AHandle, ElConstGrid* gridHandle )
{ 
    try { *gridHandle = Reinterpret(&Reinterpret(AHandle)->Grid()); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::ColConstrained() const
// ----------------------------------------------
ElError ElDistMatrixColConstrained_s
( ElConstDistMatrix_s AHandle, bool* colConst )
{
    try { *colConst = Reinterpret(AHandle)->ColConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColConstrained_d
( ElConstDistMatrix_d AHandle, bool* colConst )
{
    try { *colConst = Reinterpret(AHandle)->ColConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColConstrained_c
( ElConstDistMatrix_c AHandle, bool* colConst )
{
    try { *colConst = Reinterpret(AHandle)->ColConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColConstrained_z
( ElConstDistMatrix_z AHandle, bool* colConst )
{
    try { *colConst = Reinterpret(AHandle)->ColConstrained(); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::RowConstrained() const
// ----------------------------------------------
ElError ElDistMatrixRowConstrained_s
( ElConstDistMatrix_s AHandle, bool* rowConst )
{
    try { *rowConst = Reinterpret(AHandle)->RowConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowConstrained_d
( ElConstDistMatrix_d AHandle, bool* rowConst )
{
    try { *rowConst = Reinterpret(AHandle)->RowConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowConstrained_c
( ElConstDistMatrix_c AHandle, bool* rowConst )
{
    try { *rowConst = Reinterpret(AHandle)->RowConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowConstrained_z
( ElConstDistMatrix_z AHandle, bool* rowConst )
{
    try { *rowConst = Reinterpret(AHandle)->RowConstrained(); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::RootConstrained() const
// -----------------------------------------------
ElError ElDistMatrixRootConstrained_s
( ElConstDistMatrix_s AHandle, bool* rootConst )
{
    try { *rootConst = Reinterpret(AHandle)->RootConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRootConstrained_d
( ElConstDistMatrix_d AHandle, bool* rootConst )
{
    try { *rootConst = Reinterpret(AHandle)->RootConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRootConstrained_c
( ElConstDistMatrix_c AHandle, bool* rootConst )
{
    try { *rootConst = Reinterpret(AHandle)->RootConstrained(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRootConstrained_z
( ElConstDistMatrix_z AHandle, bool* rootConst )
{
    try { *rootConst = Reinterpret(AHandle)->RootConstrained(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::ColAlign() const
// ---------------------------------------
ElError ElDistMatrixColAlign_s( ElConstDistMatrix_s AHandle, ElInt* colAlign )
{
    try { *colAlign = Reinterpret(AHandle)->ColAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColAlign_d( ElConstDistMatrix_d AHandle, ElInt* colAlign )
{
    try { *colAlign = Reinterpret(AHandle)->ColAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColAlign_c( ElConstDistMatrix_c AHandle, ElInt* colAlign )
{
    try { *colAlign = Reinterpret(AHandle)->ColAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColAlign_z( ElConstDistMatrix_z AHandle, ElInt* colAlign )
{
    try { *colAlign = Reinterpret(AHandle)->ColAlign(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RowAlign() const
// ---------------------------------------
ElError ElDistMatrixRowAlign_s( ElConstDistMatrix_s AHandle, ElInt* rowAlign )
{
    try { *rowAlign = Reinterpret(AHandle)->RowAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowAlign_d( ElConstDistMatrix_d AHandle, ElInt* rowAlign )
{
    try { *rowAlign = Reinterpret(AHandle)->RowAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowAlign_c( ElConstDistMatrix_c AHandle, ElInt* rowAlign )
{
    try { *rowAlign = Reinterpret(AHandle)->RowAlign(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowAlign_z( ElConstDistMatrix_z AHandle, ElInt* rowAlign )
{
    try { *rowAlign = Reinterpret(AHandle)->RowAlign(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::ColShift() const
// ---------------------------------------
ElError ElDistMatrixColShift_s( ElConstDistMatrix_s AHandle, ElInt* colShift )
{
    try { *colShift = Reinterpret(AHandle)->ColShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColShift_d( ElConstDistMatrix_d AHandle, ElInt* colShift )
{
    try { *colShift = Reinterpret(AHandle)->ColShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColShift_c( ElConstDistMatrix_c AHandle, ElInt* colShift )
{
    try { *colShift = Reinterpret(AHandle)->ColShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColShift_z( ElConstDistMatrix_z AHandle, ElInt* colShift )
{
    try { *colShift = Reinterpret(AHandle)->ColShift(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RowShift() const
// ---------------------------------------
ElError ElDistMatrixRowShift_s( ElConstDistMatrix_s AHandle, ElInt* rowShift )
{
    try { *rowShift = Reinterpret(AHandle)->RowShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowShift_d( ElConstDistMatrix_d AHandle, ElInt* rowShift )
{
    try { *rowShift = Reinterpret(AHandle)->RowShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowShift_c( ElConstDistMatrix_c AHandle, ElInt* rowShift )
{
    try { *rowShift = Reinterpret(AHandle)->RowShift(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowShift_z( ElConstDistMatrix_z AHandle, ElInt* rowShift )
{
    try { *rowShift = Reinterpret(AHandle)->RowShift(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::ColRank() const
// ---------------------------------------
ElError ElDistMatrixColRank_s( ElConstDistMatrix_s AHandle, ElInt* colRank )
{
    try { *colRank = Reinterpret(AHandle)->ColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColRank_d( ElConstDistMatrix_d AHandle, ElInt* colRank )
{
    try { *colRank = Reinterpret(AHandle)->ColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColRank_c( ElConstDistMatrix_c AHandle, ElInt* colRank )
{
    try { *colRank = Reinterpret(AHandle)->ColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColRank_z( ElConstDistMatrix_z AHandle, ElInt* colRank )
{
    try { *colRank = Reinterpret(AHandle)->ColRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RowRank() const
// ---------------------------------------
ElError ElDistMatrixRowRank_s( ElConstDistMatrix_s AHandle, ElInt* rowRank )
{
    try { *rowRank = Reinterpret(AHandle)->RowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowRank_d( ElConstDistMatrix_d AHandle, ElInt* rowRank )
{
    try { *rowRank = Reinterpret(AHandle)->RowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowRank_c( ElConstDistMatrix_c AHandle, ElInt* rowRank )
{
    try { *rowRank = Reinterpret(AHandle)->RowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowRank_z( ElConstDistMatrix_z AHandle, ElInt* rowRank )
{
    try { *rowRank = Reinterpret(AHandle)->RowRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialColRank() const
// ---------------------------------------------
ElError ElDistMatrixPartialColRank_s( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColRank_d( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColRank_c( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColRank_z( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialColRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialRowRank() const
// ---------------------------------------------
ElError ElDistMatrixPartialRowRank_s( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowRank_d( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowRank_c( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowRank_z( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialRowRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialUnionColRank() const
// --------------------------------------------------
ElError ElDistMatrixPartialUnionColRank_s
( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColRank_d
( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColRank_c
( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionColRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColRank_z
( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionColRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialUnionRowRank() const
// --------------------------------------------------
ElError ElDistMatrixPartialUnionRowRank_s
( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowRank_d
( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowRank_c
( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionRowRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowRank_z
( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->PartialUnionRowRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::DistRank() const
// ---------------------------------------
ElError ElDistMatrixDistRank_s( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->DistRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistRank_d( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->DistRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistRank_c( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->DistRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistRank_z( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->DistRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::CrossRank() const
// ----------------------------------------
ElError ElDistMatrixCrossRank_s( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->CrossRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossRank_d( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->CrossRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossRank_c( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->CrossRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossRank_z( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->CrossRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RedundantRank() const
// --------------------------------------------
ElError ElDistMatrixRedundantRank_s( ElConstDistMatrix_s AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->RedundantRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantRank_d( ElConstDistMatrix_d AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->RedundantRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantRank_c( ElConstDistMatrix_c AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->RedundantRank(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantRank_z( ElConstDistMatrix_z AHandle, ElInt* rank )
{
    try { *rank = Reinterpret(AHandle)->RedundantRank(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::Root() const
// -----------------------------------
ElError ElDistMatrixRoot_s( ElConstDistMatrix_s AHandle, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->Root(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRoot_d( ElConstDistMatrix_d AHandle, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->Root(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRoot_c( ElConstDistMatrix_c AHandle, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->Root(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRoot_z( ElConstDistMatrix_z AHandle, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->Root(); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::Participating() const
// ---------------------------------------------
ElError ElDistMatrixParticipating_s
( ElConstDistMatrix_s AHandle, bool* participating )
{
    try { *participating = Reinterpret(AHandle)->Participating(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixParticipating_d
( ElConstDistMatrix_d AHandle, bool* participating )
{
    try { *participating = Reinterpret(AHandle)->Participating(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixParticipating_c
( ElConstDistMatrix_c AHandle, bool* participating )
{
    try { *participating = Reinterpret(AHandle)->Participating(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixParticipating_z
( ElConstDistMatrix_z AHandle, bool* participating )
{
    try { *participating = Reinterpret(AHandle)->Participating(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RowOwner( Int i ) const
// ----------------------------------------------
ElError ElDistMatrixRowOwner_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt* rowOwner )
{
    try { *rowOwner = Reinterpret(AHandle)->RowOwner(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowOwner_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt* rowOwner )
{
    try { *rowOwner = Reinterpret(AHandle)->RowOwner(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowOwner_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt* rowOwner )
{
    try { *rowOwner = Reinterpret(AHandle)->RowOwner(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowOwner_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt* rowOwner )
{
    try { *rowOwner = Reinterpret(AHandle)->RowOwner(i); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::ColOwner( Int j ) const
// ----------------------------------------------
ElError ElDistMatrixColOwner_s
( ElConstDistMatrix_s AHandle, ElInt j, ElInt* colOwner )
{
    try { *colOwner = Reinterpret(AHandle)->ColOwner(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColOwner_d
( ElConstDistMatrix_d AHandle, ElInt j, ElInt* colOwner )
{
    try { *colOwner = Reinterpret(AHandle)->ColOwner(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColOwner_c
( ElConstDistMatrix_c AHandle, ElInt j, ElInt* colOwner )
{
    try { *colOwner = Reinterpret(AHandle)->ColOwner(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColOwner_z
( ElConstDistMatrix_z AHandle, ElInt j, ElInt* colOwner )
{
    try { *colOwner = Reinterpret(AHandle)->ColOwner(j); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::Owner( Int i, Int j ) const
// --------------------------------------------------
ElError ElDistMatrixOwner_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt j, ElInt* owner )
{
    try { *owner = Reinterpret(AHandle)->Owner(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixOwner_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt j, ElInt* owner )
{
    try { *owner = Reinterpret(AHandle)->Owner(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixOwner_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt j, ElInt* owner )
{
    try { *owner = Reinterpret(AHandle)->Owner(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixOwner_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt j, ElInt* owner )
{
    try { *owner = Reinterpret(AHandle)->Owner(i,j); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalRow( Int i ) const
// ----------------------------------------------
ElError ElDistMatrixLocalRow_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRow_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRow_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRow_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalCol( Int j ) const
// ----------------------------------------------
ElError ElDistMatrixLocalCol_s
( ElConstDistMatrix_s AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalCol_d
( ElConstDistMatrix_d AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalCol_c
( ElConstDistMatrix_c AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalCol_z
( ElConstDistMatrix_z AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalRowOffset( Int i ) const
// ----------------------------------------------------
ElError ElDistMatrixLocalRowOffset_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRowOffset(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRowOffset_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRowOffset(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRowOffset_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRowOffset(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalRowOffset_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt* iLoc )
{
    try { *iLoc = Reinterpret(AHandle)->LocalRowOffset(i); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::LocalColOffset( Int j ) const
// ----------------------------------------------------
ElError ElDistMatrixLocalColOffset_s
( ElConstDistMatrix_s AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalColOffset(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalColOffset_d
( ElConstDistMatrix_d AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalColOffset(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalColOffset_c
( ElConstDistMatrix_c AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalColOffset(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixLocalColOffset_z
( ElConstDistMatrix_z AHandle, ElInt j, ElInt* jLoc )
{
    try { *jLoc = Reinterpret(AHandle)->LocalColOffset(j); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::GlobalRow( Int iLoc ) const
// --------------------------------------------------
ElError ElDistMatrixGlobalRow_s
( ElConstDistMatrix_s AHandle, ElInt iLoc, ElInt* i )
{
    try { *i = Reinterpret(AHandle)->GlobalRow(iLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalRow_d
( ElConstDistMatrix_d AHandle, ElInt iLoc, ElInt* i )
{
    try { *i = Reinterpret(AHandle)->GlobalRow(iLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalRow_c
( ElConstDistMatrix_c AHandle, ElInt iLoc, ElInt* i )
{
    try { *i = Reinterpret(AHandle)->GlobalRow(iLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalRow_z
( ElConstDistMatrix_z AHandle, ElInt iLoc, ElInt* i )
{
    try { *i = Reinterpret(AHandle)->GlobalRow(iLoc); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::GlobalCol( Int jLoc ) const
// --------------------------------------------------
ElError ElDistMatrixGlobalCol_s
( ElConstDistMatrix_s AHandle, ElInt jLoc, ElInt* j )
{
    try { *j = Reinterpret(AHandle)->GlobalCol(jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalCol_d
( ElConstDistMatrix_d AHandle, ElInt jLoc, ElInt* j )
{
    try { *j = Reinterpret(AHandle)->GlobalCol(jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalCol_c
( ElConstDistMatrix_c AHandle, ElInt jLoc, ElInt* j )
{
    try { *j = Reinterpret(AHandle)->GlobalCol(jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGlobalCol_z
( ElConstDistMatrix_z AHandle, ElInt jLoc, ElInt* j )
{
    try { *j = Reinterpret(AHandle)->GlobalCol(jLoc); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::IsLocalRow( Int i ) const
// -------------------------------------------------
ElError ElDistMatrixIsLocalRow_s
( ElConstDistMatrix_s AHandle, ElInt i, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalRow_d
( ElConstDistMatrix_d AHandle, ElInt i, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalRow_c
( ElConstDistMatrix_c AHandle, ElInt i, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalRow_z
( ElConstDistMatrix_z AHandle, ElInt i, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalRow(i); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::IsLocalCol( Int j ) const
// -------------------------------------------------
ElError ElDistMatrixIsLocalCol_s
( ElConstDistMatrix_s AHandle, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalCol_d
( ElConstDistMatrix_d AHandle, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalCol_c
( ElConstDistMatrix_c AHandle, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocalCol_z
( ElConstDistMatrix_z AHandle, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocalCol(j); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::IsLocal( Int i, Int j ) const
// -----------------------------------------------------
ElError ElDistMatrixIsLocal_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocal(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocal_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocal(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocal_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocal(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixIsLocal_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt j, bool* isLocal )
{
    try { *isLocal = Reinterpret(AHandle)->IsLocal(i,j); }
    CATCH
    return EL_SUCCESS;
}

// DistData DistMatrix<T,U,V>::DistData() const
// --------------------------------------------
ElError ElDistMatrixDistData_s
( ElConstDistMatrix_s AHandle, ElDistData* distData )
{
    try 
    { 
        DistData data = Reinterpret(AHandle)->DistData(); 
        *distData = Reinterpret( data );
    }
    CATCH
    return EL_SUCCESS; 
}

ElError ElDistMatrixDistData_d
( ElConstDistMatrix_d AHandle, ElDistData* distData )
{
    try 
    { 
        DistData data = Reinterpret(AHandle)->DistData(); 
        *distData = Reinterpret( data );
    }
    CATCH
    return EL_SUCCESS; 
}

ElError ElDistMatrixDistData_c
( ElConstDistMatrix_c AHandle, ElDistData* distData )
{
    try 
    { 
        DistData data = Reinterpret(AHandle)->DistData(); 
        *distData = Reinterpret( data );
    }
    CATCH
    return EL_SUCCESS; 
}

ElError ElDistMatrixDistData_z
( ElConstDistMatrix_z AHandle, ElDistData* distData )
{
    try 
    { 
        DistData data = Reinterpret(AHandle)->DistData(); 
        *distData = Reinterpret( data );
    }
    CATCH
    return EL_SUCCESS; 
}

// mpi::Comm DistMatrix<T,U,V>::DistComm() const
// ---------------------------------------------
ElError ElDistMatrixDistComm_s( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->DistComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistComm_d( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->DistComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistComm_c( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->DistComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistComm_z( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->DistComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::CrossComm() const
// ----------------------------------------------
ElError ElDistMatrixCrossComm_s( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->CrossComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossComm_d( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->CrossComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossComm_c( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->CrossComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossComm_z( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->CrossComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::RedundantComm() const
// --------------------------------------------------
ElError ElDistMatrixRedundantComm_s
( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RedundantComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantComm_d
( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RedundantComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantComm_c
( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RedundantComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantComm_z
( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RedundantComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::ColComm() const
// --------------------------------------------
ElError ElDistMatrixColComm_s( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->ColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColComm_d( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->ColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColComm_c( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->ColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColComm_z( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->ColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::RowComm() const
// --------------------------------------------
ElError ElDistMatrixRowComm_s( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowComm_d( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowComm_c( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowComm_z( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->RowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::PartialColComm() const
// ---------------------------------------------------
ElError ElDistMatrixPartialColComm_s
( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColComm_d
( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColComm_c
( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColComm_z
( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::PartialRowComm() const
// ---------------------------------------------------
ElError ElDistMatrixPartialRowComm_s
( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowComm_d
( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowComm_c
( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowComm_z
( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::PartialUnionColComm() const
// --------------------------------------------------------
ElError ElDistMatrixPartialUnionColComm_s
( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColComm_d
( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColComm_c
( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColComm_z
( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionColComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// mpi::Comm DistMatrix<T,U,V>::PartialUnionRowComm() const
// --------------------------------------------------------
ElError ElDistMatrixPartialUnionRowComm_s
( ElConstDistMatrix_s AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowComm_d
( ElConstDistMatrix_d AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowComm_c
( ElConstDistMatrix_c AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowComm_z
( ElConstDistMatrix_z AHandle, MPI_Comm* comm )
{
    try { *comm = Reinterpret(AHandle)->PartialUnionRowComm().comm; }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::ColStride() const
// ----------------------------------------
ElError ElDistMatrixColStride_s( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->ColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColStride_d( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->ColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColStride_c( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->ColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixColStride_z( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->ColStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RowStride() const
// ----------------------------------------
ElError ElDistMatrixRowStride_s( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->RowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowStride_d( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->RowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowStride_c( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->RowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRowStride_z( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->RowStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialColStride() const
// -----------------------------------------------
ElError ElDistMatrixPartialColStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialColStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialColStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialRowStride() const
// -----------------------------------------------
ElError ElDistMatrixPartialRowStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialRowStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialRowStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialUnionColStride() const
// ----------------------------------------------------
ElError ElDistMatrixPartialUnionColStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionColStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionColStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionColStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::PartialUnionRowStride() const
// ----------------------------------------------------
ElError ElDistMatrixPartialUnionRowStride_s
( ElConstDistMatrix_s AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowStride_d
( ElConstDistMatrix_d AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowStride_c
( ElConstDistMatrix_c AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionRowStride(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixPartialUnionRowStride_z
( ElConstDistMatrix_z AHandle, ElInt* stride )
{
    try { *stride = Reinterpret(AHandle)->PartialUnionRowStride(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::DistSize() const
// ---------------------------------------
ElError ElDistMatrixDistSize_s( ElConstDistMatrix_s AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->DistSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistSize_d( ElConstDistMatrix_d AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->DistSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistSize_c( ElConstDistMatrix_c AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->DistSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDistSize_z( ElConstDistMatrix_z AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->DistSize(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::CrossSize() const
// ----------------------------------------
ElError ElDistMatrixCrossSize_s( ElConstDistMatrix_s AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->CrossSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossSize_d( ElConstDistMatrix_d AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->CrossSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossSize_c( ElConstDistMatrix_c AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->CrossSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixCrossSize_z( ElConstDistMatrix_z AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->CrossSize(); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::RedundantSize() const
// --------------------------------------------
ElError ElDistMatrixRedundantSize_s
( ElConstDistMatrix_s AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->RedundantSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantSize_d
( ElConstDistMatrix_d AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->RedundantSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantSize_c
( ElConstDistMatrix_c AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->RedundantSize(); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixRedundantSize_z
( ElConstDistMatrix_z AHandle, ElInt* commSize )
{
    try { *commSize = Reinterpret(AHandle)->RedundantSize(); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Get( Int i, Int j ) const
// -------------------------------------------------
ElError ElDistMatrixGet_s
( ElConstDistMatrix_s AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = Reinterpret(AHandle)->Get(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGet_d
( ElConstDistMatrix_d AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = Reinterpret(AHandle)->Get(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGet_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt j, complex_float* val )
{
    try 
    { 
        Complex<float> alpha = Reinterpret(AHandle)->Get(i,j);
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGet_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt j, complex_double* val )
{
    try 
    { 
        Complex<double> alpha = Reinterpret(AHandle)->Get(i,j); 
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::GetRealPart( Int i, Int j ) const
// ---------------------------------------------------------
ElError ElDistMatrixGetRealPart_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = Reinterpret(AHandle)->GetRealPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetRealPart_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = Reinterpret(AHandle)->GetRealPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::GetImagPart( Int i, Int j ) const
// ---------------------------------------------------------
ElError ElDistMatrixGetImagPart_c
( ElConstDistMatrix_c AHandle, ElInt i, ElInt j, float* val )
{
    try { *val = Reinterpret(AHandle)->GetImagPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetImagPart_z
( ElConstDistMatrix_z AHandle, ElInt i, ElInt j, double* val )
{
    try { *val = Reinterpret(AHandle)->GetImagPart(i,j); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Set( Int i, Int j, T alpha )
// ----------------------------------------------------
ElError ElDistMatrixSet_s
( ElDistMatrix_s AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->Set(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSet_d
( ElDistMatrix_d AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->Set(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSet_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { Reinterpret(AHandle)->Set(i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSet_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { Reinterpret(AHandle)->Set(i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetRealPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
ElError ElDistMatrixSetRealPart_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRealPart_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->SetRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetImagPart( Int i, Int j, Base<T> alpha )
// ------------------------------------------------------------------
ElError ElDistMatrixSetImagPart_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetImagPart_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->SetImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Update( Int i, Int j, T alpha )
// -------------------------------------------------------
ElError ElDistMatrixUpdate_s
( ElDistMatrix_s AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->Update(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdate_d
( ElDistMatrix_d AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->Update(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdate_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, complex_float alpha )
{
    try { Reinterpret(AHandle)->Update
          (i,j,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdate_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, complex_double alpha )
{
    try { Reinterpret(AHandle)->Update
          (i,j,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateRealPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
ElError ElDistMatrixUpdateRealPart_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateRealPart_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->UpdateRealPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateImagPart( Int i, Int j, Base<T> alpha )
// ---------------------------------------------------------------------
ElError ElDistMatrixUpdateImagPart_c
( ElDistMatrix_c AHandle, ElInt i, ElInt j, float alpha )
{
    try { Reinterpret(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateImagPart_z
( ElDistMatrix_z AHandle, ElInt i, ElInt j, double alpha )
{
    try { Reinterpret(AHandle)->UpdateImagPart(i,j,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeReal( Int i, Int j )
// ------------------------------------------------
ElError ElDistMatrixMakeReal_c( ElDistMatrix_c AHandle, ElInt i, ElInt j )
{
    try { Reinterpret(AHandle)->MakeReal(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeReal_z( ElDistMatrix_z AHandle, ElInt i, ElInt j )
{
    try { Reinterpret(AHandle)->MakeReal(i,j); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::Conjugate( Int i, Int j )
// -------------------------------------------------
ElError ElDistMatrixConjugate_c( ElDistMatrix_c AHandle, ElInt i, ElInt j )
{
    try { Reinterpret(AHandle)->Conjugate(i,j); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixConjugate_z( ElDistMatrix_z AHandle, ElInt i, ElInt j )
{
    try { Reinterpret(AHandle)->Conjugate(i,j); }
    CATCH
    return EL_SUCCESS;
}

// T DistMatrix<T,U,V>::GetLocal( Int iLoc, Int jLoc ) const
// ---------------------------------------------------------
ElError ElDistMatrixGetLocal_s
( ElConstDistMatrix_s AHandle, ElInt iLoc, ElInt jLoc, float* val )
{
    try { *val = Reinterpret(AHandle)->GetLocal(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocal_d
( ElConstDistMatrix_d AHandle, ElInt iLoc, ElInt jLoc, double* val )
{
    try { *val = Reinterpret(AHandle)->GetLocal(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocal_c
( ElConstDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, complex_float* val )
{
    try 
    { 
        Complex<float> alpha = Reinterpret(AHandle)->GetLocal(iLoc,jLoc);
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocal_z
( ElConstDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, complex_double* val )
{
    try 
    { 
        Complex<double> alpha = Reinterpret(AHandle)->GetLocal(iLoc,jLoc);
        val->real = alpha.real();
        val->imag = alpha.imag();
    }
    CATCH
    return EL_SUCCESS;
}

// Base<T> DistMatrix<T,U,V>::GetLocalRealPart( Int iLoc, Int jLoc ) const
// -----------------------------------------------------------------------
ElError ElDistMatrixGetLocalRealPart_c
( ElConstDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float* val )
{
    try { *val = Reinterpret(AHandle)->GetLocalRealPart(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocalRealPart_z
( ElConstDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double* val )
{
    try { *val = Reinterpret(AHandle)->GetLocalRealPart(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

// Base<T> DistMatrix<T,U,V>::GetLocalImagPart( Int iLoc, Int jLoc ) const
// -----------------------------------------------------------------------
ElError ElDistMatrixGetLocalImagPart_c
( ElConstDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float* val )
{
    try { *val = Reinterpret(AHandle)->GetLocalImagPart(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocalImagPart_z
( ElConstDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double* val )
{
    try { *val = Reinterpret(AHandle)->GetLocalImagPart(iLoc,jLoc); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetLocal( Int iLoc, Int jLoc, T alpha )
// ---------------------------------------------------------------
ElError ElDistMatrixSetLocal_s
( ElDistMatrix_s AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->SetLocal(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocal_d
( ElDistMatrix_d AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->SetLocal(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocal_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, complex_float alpha )
{
    try { Reinterpret(AHandle)->SetLocal
          (iLoc,jLoc,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocal_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, complex_double alpha )
{
    try { Reinterpret(AHandle)->SetLocal
          (iLoc,jLoc,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha )
// -----------------------------------------------------------------------------
ElError ElDistMatrixSetLocalRealPart_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->SetLocalRealPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocalRealPart_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->SetLocalRealPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha )
// -----------------------------------------------------------------------------
ElError ElDistMatrixSetLocalImagPart_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->SetLocalImagPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocalImagPart_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->SetLocalImagPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateLocal( Int iLoc, Int jLoc, T alpha )
// ------------------------------------------------------------------
ElError ElDistMatrixUpdateLocal_s
( ElDistMatrix_s AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->UpdateLocal(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocal_d
( ElDistMatrix_d AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->UpdateLocal(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocal_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, complex_float alpha )
{
    try { Reinterpret(AHandle)->UpdateLocal
          (iLoc,jLoc,Complex<float>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocal_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, complex_double alpha )
{
    try { Reinterpret(AHandle)->UpdateLocal
          (iLoc,jLoc,Complex<double>(alpha.real,alpha.imag)); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateLocalRealPart
// ( Int iLoc, Int jLoc, Base<T> alpha )
// -------------------------------------------
ElError ElDistMatrixUpdateLocalRealPart_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->UpdateLocalRealPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocalRealPart_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->UpdateLocalRealPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateLocalImagPart
// ( Int iLoc, Int jLoc, Base<T> alpha )
// -------------------------------------------
ElError ElDistMatrixUpdateLocalImagPart_c
( ElDistMatrix_c AHandle, ElInt iLoc, ElInt jLoc, float alpha )
{
    try { Reinterpret(AHandle)->UpdateLocalImagPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocalImagPart_z
( ElDistMatrix_z AHandle, ElInt iLoc, ElInt jLoc, double alpha )
{
    try { Reinterpret(AHandle)->UpdateLocalImagPart(iLoc,jLoc,alpha); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeDiagonalReal( Int offset )
// ------------------------------------------------------
ElError ElDistMatrixMakeDiagonalReal_c( ElDistMatrix_c AHandle, ElInt offset )
{
    try { Reinterpret(AHandle)->MakeDiagonalReal(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeDiagonalReal_z( ElDistMatrix_z AHandle, ElInt offset )
{
    try { Reinterpret(AHandle)->MakeDiagonalReal(offset); }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::ConjugateDiagonal( Int offset )
// -------------------------------------------------------
ElError ElDistMatrixConjugateDiagonal_c( ElDistMatrix_c AHandle, ElInt offset )
{
    try { Reinterpret(AHandle)->ConjugateDiagonal(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixConjugateDiagonal_z( ElDistMatrix_z AHandle, ElInt offset )
{
    try { Reinterpret(AHandle)->ConjugateDiagonal(offset); }
    CATCH
    return EL_SUCCESS;
}

// bool DistMatrix<T,U,V>::DiagonalAlignedWith
// ( const DistData& data, Int offset ) const
// -------------------------------------------
ElError ElDistMatrixDiagonalAlignedWith_s
( ElConstDistMatrix_s AHandle, ElDistData distData, ElInt offset, 
  bool* aligned )
{
    try { *aligned = Reinterpret(AHandle)->DiagonalAlignedWith
                     (Reinterpret(distData),offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlignedWith_d
( ElConstDistMatrix_d AHandle, ElDistData distData, ElInt offset, 
  bool* aligned )
{
    try { *aligned = Reinterpret(AHandle)->DiagonalAlignedWith
                     (Reinterpret(distData),offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlignedWith_c
( ElConstDistMatrix_c AHandle, ElDistData distData, ElInt offset, 
  bool* aligned )
{
    try { *aligned = Reinterpret(AHandle)->DiagonalAlignedWith
                     (Reinterpret(distData),offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlignedWith_z
( ElConstDistMatrix_z AHandle, ElDistData distData, ElInt offset, 
  bool* aligned )
{
    try { *aligned = Reinterpret(AHandle)->DiagonalAlignedWith
                     (Reinterpret(distData),offset); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::DiagonalRoot( Int offset ) const
// -------------------------------------------------------
ElError ElDistMatrixDiagonalRoot_s
( ElConstDistMatrix_s AHandle, ElInt offset, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->DiagonalRoot(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalRoot_d
( ElConstDistMatrix_d AHandle, ElInt offset, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->DiagonalRoot(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalRoot_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->DiagonalRoot(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalRoot_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElInt* root )
{
    try { *root = Reinterpret(AHandle)->DiagonalRoot(offset); }
    CATCH
    return EL_SUCCESS;
}

// Int DistMatrix<T,U,V>::DiagonalAlign( Int offset ) const
// --------------------------------------------------------
ElError ElDistMatrixDiagonalAlign_s
( ElConstDistMatrix_s AHandle, ElInt offset, ElInt* align )
{
    try { *align = Reinterpret(AHandle)->DiagonalAlign(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlign_d
( ElConstDistMatrix_d AHandle, ElInt offset, ElInt* align )
{
    try { *align = Reinterpret(AHandle)->DiagonalAlign(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlign_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElInt* align )
{
    try { *align = Reinterpret(AHandle)->DiagonalAlign(offset); }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixDiagonalAlign_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElInt* align )
{
    try { *align = Reinterpret(AHandle)->DiagonalAlign(offset); }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"

// DistMatrix<T,UDiag,VDiag> DistMatrix<T,U,V>::GetDiagonal( Int offset ) const
// ----------------------------------------------------------------------------
template<typename T,Dist CDist,Dist RDist>
void ElDistMatrixGetDiagonalKernel
( const AbstractDistMatrix<T>* AAbs, Int offset, AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    if( U == CDist && V == RDist )
    {
        auto A = dynamic_cast<const DistMatrix<T,CDist,RDist>*>(AAbs);
        DynamicCastCheck(A);
        auto* d = new DistMatrix<T,DiagColDist<CDist,RDist>(),
                                   DiagRowDist<CDist,RDist>()>(g);
        A->GetDiagonal( *d, offset );
        *dAbs = d;
    }
}

template<typename T>
ElError ElDistMatrixGetDiagonal
( const AbstractDistMatrix<T>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    try 
    {
        ElDistMatrixGetDiagonalKernel<T,CIRC,CIRC>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,MC,  MR  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,MC,  STAR>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,MD,  STAR>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,MR,  MC  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,MR,  STAR>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,MC  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,MD  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,MR  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,STAR>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,VC  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,STAR,VR  >( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,VC,  STAR>( AAbs, offset, dAbs );
        ElDistMatrixGetDiagonalKernel<T,VR,  STAR>( AAbs, offset, dAbs );
    }
    CATCH
    return EL_SUCCESS;
}

extern "C" { 

ElError ElDistMatrixGetDiagonal_s
( ElConstDistMatrix_s AHandle, ElInt offset, ElDistMatrix_s* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<float>* dAbs;
    ElError error = ElDistMatrixGetDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

ElError ElDistMatrixGetDiagonal_d
( ElConstDistMatrix_d AHandle, ElInt offset, ElDistMatrix_d* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<double>* dAbs;
    ElError error = ElDistMatrixGetDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

ElError ElDistMatrixGetDiagonal_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElDistMatrix_c* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<Complex<float>>* dAbs;
    ElError error = ElDistMatrixGetDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

ElError ElDistMatrixGetDiagonal_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElDistMatrix_z* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<Complex<double>>* dAbs;
    ElError error = ElDistMatrixGetDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

} // extern "C"

// DistMatrix<Base<T>,UDiag,VDiag> 
// DistMatrix<T,U,V>::GetRealPartOfDiagonal( Int offset ) const
// ------------------------------------------------------------
template<typename T,Dist CDist,Dist RDist>
void ElDistMatrixGetRealPartOfDiagonalKernel
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    if( U == CDist && V == RDist )
    {
        auto A = dynamic_cast<const DistMatrix<Complex<T>,CDist,RDist>*>(AAbs);
        DynamicCastCheck(A);
        auto* d = new DistMatrix<T,DiagColDist<CDist,RDist>(),
                                   DiagRowDist<CDist,RDist>()>(g);
        A->GetRealPartOfDiagonal( *d, offset );
        *dAbs = d;
    }
}

template<typename T>
ElError ElDistMatrixGetRealPartOfDiagonal
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    try 
    {
        ElDistMatrixGetRealPartOfDiagonalKernel<T,CIRC,CIRC>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,MC,  MR  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,MC,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,MD,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,MR,  MC  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,MR,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,MC  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,MD  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,MR  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,VC  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,STAR,VR  >(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,VC,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetRealPartOfDiagonalKernel<T,VR,  STAR>(AAbs,offset,dAbs);
    }
    CATCH
    return EL_SUCCESS;
}

extern "C" {

ElError ElDistMatrixGetRealPartOfDiagonal_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElDistMatrix_s* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<float>* dAbs;
    ElError error = ElDistMatrixGetRealPartOfDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

ElError ElDistMatrixGetRealPartOfDiagonal_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElDistMatrix_d* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<double>* dAbs;
    ElError error = ElDistMatrixGetRealPartOfDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

} // extern "C"

// DistMatrix<Base<T>,UDiag,VDiag> 
// DistMatrix<T,U,V>::GetImagPartOfDiagonal( Int offset ) const
// ------------------------------------------------------------
template<typename T,Dist CDist,Dist RDist>
void ElDistMatrixGetImagPartOfDiagonalKernel
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    const Grid& g = AAbs->Grid();
    Dist U = AAbs->DistData().colDist;
    Dist V = AAbs->DistData().rowDist;
    if( U == CDist && V == RDist )
    {
        auto A = dynamic_cast<const DistMatrix<Complex<T>,CDist,RDist>*>(AAbs);
        DynamicCastCheck(A);
        auto* d = new DistMatrix<T,DiagColDist<CDist,RDist>(),
                                   DiagRowDist<CDist,RDist>()>(g);
        A->GetImagPartOfDiagonal( *d, offset );
        *dAbs = d;
    }
}

template<typename T>
ElError ElDistMatrixGetImagPartOfDiagonal
( const AbstractDistMatrix<Complex<T>>* AAbs, Int offset, 
        AbstractDistMatrix<T>** dAbs )
{
    try 
    {
        ElDistMatrixGetImagPartOfDiagonalKernel<T,CIRC,CIRC>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,MC,  MR  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,MC,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,MD,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,MR,  MC  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,MR,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,MC  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,MD  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,MR  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,VC  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,STAR,VR  >(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,VC,  STAR>(AAbs,offset,dAbs);
        ElDistMatrixGetImagPartOfDiagonalKernel<T,VR,  STAR>(AAbs,offset,dAbs);
    }
    CATCH
    return EL_SUCCESS;
}

extern "C" {

ElError ElDistMatrixGetImagPartOfDiagonal_c
( ElConstDistMatrix_c AHandle, ElInt offset, ElDistMatrix_s* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<float>* dAbs;
    ElError error = ElDistMatrixGetImagPartOfDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

ElError ElDistMatrixGetImagPartOfDiagonal_z
( ElConstDistMatrix_z AHandle, ElInt offset, ElDistMatrix_d* dHandle )
{
    auto AAbs = Reinterpret(AHandle);
    AbstractDistMatrix<double>* dAbs;
    ElError error = ElDistMatrixGetImagPartOfDiagonal( AAbs, offset, &dAbs );
    *dHandle = Reinterpret(dAbs);
    return error;
}

// TODO: More diagonal manipulation
// ================================

// DistMatrix<T,STAR,STAR> DistMatrix<T,U,V>::GetSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElDistMatrixGetSubmatrix_s
( ElConstDistMatrix_s AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<float,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetSubmatrix_d
( ElConstDistMatrix_d AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<double,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_c* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<Complex<float>,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_c)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_z* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<Complex<double>,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_z)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetRealPartOfSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElDistMatrixGetRealPartOfSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<float,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetRealPartOfSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<double,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// DistMatrix<Base<T>,STAR,STAR> DistMatrix<T,U,V>::GetImagPartOfSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElDistMatrixGetImagPartOfSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<float,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetImagPartOfSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowInds,
  ElInt numColInds, const ElInt* colInds, ElDistMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new DistMatrix<double,STAR,STAR>(A->Grid());

        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->GetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElDistMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   const DistMatrix<T,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixSetSubmatrix_s
( ElDistMatrix_s AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetSubmatrix_d
( ElDistMatrix_d AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_c ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = 
            dynamic_cast<const DistMatrix<Complex<float>,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_z ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = 
            dynamic_cast<const DistMatrix<Complex<double>,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetRealPartOfSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   const DistMatrix<Base<T>,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixSetRealPartOfSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRealPartOfSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetRealPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetImagPartOfSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   const DistMatrix<Base<T>,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixSetImagPartOfSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetImagPartOfSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->SetImagPartOfSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   T alpha, const DistMatrix<T,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixUpdateSubmatrix_s
( ElDistMatrix_s AHandle, const ElInt* rowInds, const ElInt* colInds,
  float alpha, ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateSubmatrix_d
( ElDistMatrix_d AHandle, const ElInt* rowInds, const ElInt* colInds,
  double alpha, ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  complex_float alpha, ElConstDistMatrix_c ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = 
            dynamic_cast<const DistMatrix<Complex<float>,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateSubmatrix
        ( rowIndVec, colIndVec, Complex<float>(alpha.real,alpha.imag), *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  complex_double alpha, ElConstDistMatrix_z ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = 
            dynamic_cast<const DistMatrix<Complex<double>,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateSubmatrix
        ( rowIndVec, colIndVec, Complex<double>(alpha.real,alpha.imag), *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateRealPartOfSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixUpdateRealPartOfSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  float alpha, ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateRealPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateRealPartOfSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  double alpha, ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateRealPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateImagPartOfSubmatrix
// ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
//   Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub );
// -----------------------------------------------------------------
ElError ElDistMatrixUpdateImagPartOfSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowInds, const ElInt* colInds,
  float alpha, ElConstDistMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<float,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateImagPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateImagPartOfSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowInds, const ElInt* colInds,
  double alpha, ElConstDistMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASubADM = Reinterpret(ASubHandle);
        auto ASub = dynamic_cast<const DistMatrix<double,STAR,STAR>*>(ASubADM);
        DynamicCastCheck(ASub); 

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        A->UpdateImagPartOfSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeSubmatrixReal
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
// --------------------------------------------------------------------
ElError ElDistMatrixMakeSubmatrixReal_c
( ElDistMatrix_c AHandle, ElInt numRowInds, const ElInt* rowInds,
                          ElInt numColInds, const ElInt* colInds )
{
    try 
    {
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        Reinterpret(AHandle)->MakeSubmatrixReal( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeSubmatrixReal_z
( ElDistMatrix_z AHandle, ElInt numRowInds, const ElInt* rowInds,
                          ElInt numColInds, const ElInt* colInds )
{
    try 
    {
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        Reinterpret(AHandle)->MakeSubmatrixReal( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::ConjugateSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds )
// --------------------------------------------------------------------
ElError ElDistMatrixConjugateSubmatrix_c
( ElDistMatrix_c AHandle, ElInt numRowInds, const ElInt* rowInds,
                          ElInt numColInds, const ElInt* colInds )
{
    try 
    {
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        Reinterpret(AHandle)->ConjugateSubmatrix( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixConjugateSubmatrix_z
( ElDistMatrix_z AHandle, ElInt numRowInds, const ElInt* rowInds,
                          ElInt numColInds, const ElInt* colInds )
{
    try 
    {
        std::vector<Int> rowIndVec(rowInds,rowInds+numRowInds),
                         colIndVec(colInds,colInds+numColInds);
        Reinterpret(AHandle)->ConjugateSubmatrix( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

// Matrix<T> DistMatrix<T,U,V>::GetLocalSubmatrix
// ( const std::vector<Int>& rowIndsLoc, 
//   const std::vector<Int>& colIndsLoc ) const
// ------------------------------------------------------------
ElError ElDistMatrixGetLocalSubmatrix_s
( ElConstDistMatrix_s AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<float>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocalSubmatrix_d
( ElConstDistMatrix_d AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<double>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocalSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_c* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<Complex<float>>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_c)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetLocalSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_z* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<Complex<double>>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_z)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// Matrix<Base<T>> DistMatrix<T,U,V>::GetRealPartOfLocalSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElDistMatrixGetRealPartOfLocalSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<float>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetRealPartOfLocalSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<double>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// Matrix<Base<T>> DistMatrix<T,U,V>::GetImagPartOfLocalSubmatrix
// ( const std::vector<Int>& rowInds, const std::vector<Int>& colInds ) const
// --------------------------------------------------------------------------
ElError ElDistMatrixGetImagPartOfLocalSubmatrix_c
( ElConstDistMatrix_c AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_s* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<float>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_s)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixGetImagPartOfLocalSubmatrix_z
( ElConstDistMatrix_z AHandle,
  ElInt numRowInds, const ElInt* rowIndsLoc,
  ElInt numColInds, const ElInt* colIndsLoc, ElMatrix_d* ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);
        auto ASub = new Matrix<double>;

        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->GetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
        *ASubHandle = (ElMatrix_d)ASub;
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<T>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixSetLocalSubmatrix_s
( ElDistMatrix_s AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocalSubmatrix_d
( ElDistMatrix_d AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_c ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_z ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetRealPartOfLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<Base<T>>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixSetRealPartOfLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetRealPartOfLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetRealPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::SetImagPartOfLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<Base<T>>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixSetImagPartOfLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixSetImagPartOfLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->SetImagPartOfLocalSubmatrix( rowIndVec, colIndVec, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<T>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixUpdateLocalSubmatrix_s
( ElDistMatrix_s AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  float alpha, ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocalSubmatrix_d
( ElDistMatrix_d AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  double alpha, ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  complex_float alpha, ElConstMatrix_c ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateLocalSubmatrix
        ( rowIndVec, colIndVec, Complex<float>(alpha.real,alpha.imag), *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  complex_double alpha, ElConstMatrix_z ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateLocalSubmatrix
        ( rowIndVec, colIndVec, Complex<double>(alpha.real,alpha.imag), *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateRealPartOfLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<Base<T>>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  float alpha, ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateRealPartOfLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateRealPartOfLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  double alpha, ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateRealPartOfLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::UpdateImagPartOfLocalSubmatrix
// ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
//   const Matrix<Base<T>>& ASub );
// -----------------------------------------------------------------------
ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_c
( ElDistMatrix_c AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  float alpha, ElConstMatrix_s ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateImagPartOfLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixUpdateImagPartOfLocalSubmatrix_z
( ElDistMatrix_z AHandle, const ElInt* rowIndsLoc, const ElInt* colIndsLoc,
  double alpha, ElConstMatrix_d ASubHandle )
{
    try
    {
        auto A = Reinterpret(AHandle);    
        auto ASub = Reinterpret(ASubHandle);

        const Int numRowInds = ASub->Height();
        const Int numColInds = ASub->Width();
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        A->UpdateImagPartOfLocalSubmatrix( rowIndVec, colIndVec, alpha, *ASub );
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::MakeLocalSubmatrixReal
// ( const std::vector<Int>& rowIndsLoc, const std::vector<Int>& colIndsLoc )
// --------------------------------------------------------------------------
ElError ElDistMatrixMakeLocalSubmatrixReal_c
( ElDistMatrix_c AHandle, ElInt numRowInds, const ElInt* rowIndsLoc,
                          ElInt numColInds, const ElInt* colIndsLoc )
{
    try 
    {
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        Reinterpret(AHandle)->MakeLocalSubmatrixReal( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixMakeLocalSubmatrixReal_z
( ElDistMatrix_z AHandle, ElInt numRowInds, const ElInt* rowIndsLoc,
                          ElInt numColInds, const ElInt* colIndsLoc )
{
    try 
    {
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        Reinterpret(AHandle)->MakeLocalSubmatrixReal( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

// void DistMatrix<T,U,V>::ConjugateLocalSubmatrix
// ( const std::vector<Int>& rowIndsLoc, const std::vector<Int>& colIndsLoc )
// --------------------------------------------------------------------------
ElError ElDistMatrixConjugateLocalSubmatrix_c
( ElDistMatrix_c AHandle, ElInt numRowInds, const ElInt* rowIndsLoc,
                          ElInt numColInds, const ElInt* colIndsLoc )
{
    try 
    {
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        Reinterpret(AHandle)->ConjugateLocalSubmatrix( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

ElError ElDistMatrixConjugateLocalSubmatrix_z
( ElDistMatrix_z AHandle, ElInt numRowInds, const ElInt* rowIndsLoc,
                          ElInt numColInds, const ElInt* colIndsLoc )
{
    try 
    {
        std::vector<Int> rowIndVec(rowIndsLoc,rowIndsLoc+numRowInds),
                         colIndVec(colIndsLoc,colIndsLoc+numColInds);
        Reinterpret(AHandle)->ConjugateLocalSubmatrix( rowIndVec, colIndVec ); 
    }
    CATCH
    return EL_SUCCESS;
}

} // extern "C"
