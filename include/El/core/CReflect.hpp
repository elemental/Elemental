/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_CREFLECT_C_HPP
#define EL_CORE_CREFLECT_C_HPP

#define EL_CATCH \
  catch( std::bad_alloc& e ) \
  { El::ReportException(e); return EL_ALLOC_ERROR; } \
  catch( El::ArgException& e ) \
  { El::ReportException(e); return EL_ARG_ERROR; } \
  catch( std::logic_error& e ) \
  { El::ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { El::ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { El::ReportException(e); return EL_ERROR; }

#define EL_TRY(payload) \
  try { payload; } EL_CATCH \
  return EL_SUCCESS;

#define EL_RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

namespace El {

template<typename T>
struct CReflectType { typedef T type; };

// ElInt and Int are typedef's
/*
template<> struct CReflectType<ElInt> { typedef Int type; };
template<> struct CReflectType<Int> { typedef ElInt type; };
*/

template<> struct CReflectType<complex_float> 
{ typedef Complex<float> type; };
template<> struct CReflectType<complex_double> 
{ typedef Complex<double> type; };

template<> struct CReflectType<Complex<float>> 
{ typedef complex_float type; };
template<> struct CReflectType<Complex<double>> 
{ typedef complex_double type; };

#define CREFLECT(T) typename CReflectType<T>::type

template<typename T>
inline void DynamicCastCheck( T* A )
{ if( A == nullptr ) RuntimeError("Dynamic cast failed"); }

inline string CReflect( const char* name ) { return string(name); }
// NOTE: This creates a deep copy and the pointer should be deleted later
inline char* CReflect( const string& name ) 
{
    const auto size = name.size();
    char* buffer = new char[size+1];
    memcpy( buffer, name.c_str(), size+1 );
    return buffer;
}

inline Range<Int> CReflect( ElRange_i rangeC )
{ return Range<Int>(rangeC.beg,rangeC.end); }
inline ElRange_i CReflect( Range<Int> range )
{ 
    ElRange_i rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Range<float> CReflect( ElRange_s rangeC )
{ return Range<float>(rangeC.beg,rangeC.end); }
inline ElRange_s CReflect( Range<float> range )
{ 
    ElRange_s rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Range<double> CReflect( ElRange_d rangeC )
{ return Range<double>(rangeC.beg,rangeC.end); }
inline ElRange_d CReflect( Range<double> range )
{ 
    ElRange_d rangeC; 
    rangeC.beg = range.beg; 
    rangeC.end = range.end; 
    return rangeC; 
}

inline Orientation CReflect( ElOrientation orient ) 
{ return static_cast<Orientation>(orient); }
inline ElOrientation CReflect( Orientation orient )
{ return static_cast<ElOrientation>(orient); }

inline LeftOrRight CReflect( ElLeftOrRight side )
{ return static_cast<LeftOrRight>(side); }
inline ElLeftOrRight CReflect( LeftOrRight side )
{ return static_cast<ElLeftOrRight>(side); }

inline UpperOrLower CReflect( ElUpperOrLower uplo )
{ return static_cast<UpperOrLower>(uplo); }
inline ElUpperOrLower CReflect( UpperOrLower uplo )
{ return static_cast<ElUpperOrLower>(uplo); }

inline UnitOrNonUnit CReflect( ElUnitOrNonUnit diag )
{ return static_cast<UnitOrNonUnit>(diag); }
inline ElUnitOrNonUnit CReflect( UnitOrNonUnit diag )
{ return static_cast<ElUnitOrNonUnit>(diag); }

inline VerticalOrHorizontal CReflect( ElVerticalOrHorizontal dir )
{ return static_cast<VerticalOrHorizontal>(dir); }
inline ElVerticalOrHorizontal CReflect( VerticalOrHorizontal dir )
{ return static_cast<ElVerticalOrHorizontal>(dir); }

inline ForwardOrBackward CReflect( ElForwardOrBackward order )
{ return static_cast<ForwardOrBackward>(order); }
inline ElForwardOrBackward CReflect( ForwardOrBackward order )
{ return static_cast<ElForwardOrBackward>(order); }

inline Conjugation CReflect( ElConjugation conjugation )
{ return static_cast<Conjugation>(conjugation); }
inline ElConjugation CReflect( Conjugation conjugation )
{ return static_cast<ElConjugation>(conjugation); }

// Dist
// ----
inline Dist   CReflect( ElDist dist ) { return static_cast<  Dist>(dist); }
inline ElDist CReflect(   Dist dist ) { return static_cast<ElDist>(dist); }

// Grid
// ----
inline   GridOrder     CReflect( ElGridOrderType order )
{ return static_cast<  GridOrder    >(order); }
inline ElGridOrderType CReflect(   GridOrder     order )
{ return static_cast<ElGridOrderType>(order); }

inline Grid* CReflect( ElGrid grid )
{ return EL_RC(Grid*,grid); }
inline ElGrid CReflect( Grid* grid )
{ return (ElGrid)EL_RC(struct ElGrid_sDummy*,grid); }

inline const Grid* CReflect( ElConstGrid grid )
{ return EL_RC(const Grid*,grid); }
inline ElConstGrid CReflect( const Grid* grid )
{ return (ElConstGrid)EL_RC(const struct ElGrid_sDummy*,grid); }

// Complex<T>
// ----------
inline complex_float* CReflect( Complex<float>* buffer )
{ return EL_RC(complex_float*,buffer); }

inline complex_double* CReflect( Complex<double>* buffer )
{ return EL_RC(complex_double*,buffer); }

inline const complex_float* CReflect( const Complex<float>* buffer )
{ return EL_RC(const complex_float*,buffer); }

inline const complex_double* CReflect( const Complex<double>* buffer )
{ return EL_RC(const complex_double*,buffer); }

inline Complex<float>* CReflect( complex_float* buffer )
{ return EL_RC(Complex<float>*,buffer); }

inline Complex<double>* CReflect( complex_double* buffer )
{ return EL_RC(Complex<double>*,buffer); }

inline const Complex<float>* CReflect( const complex_float* buffer )
{ return EL_RC(const Complex<float>*,buffer); }

inline const Complex<double>* CReflect( const complex_double* buffer )
{ return EL_RC(const Complex<double>*,buffer); }

inline Complex<float> CReflect( complex_float alpha )
{ return Complex<float>(alpha.real,alpha.imag); }

inline Complex<double> CReflect( complex_double alpha )
{ return Complex<double>(alpha.real,alpha.imag); }

inline complex_float CReflect( Complex<float> alpha )
{ complex_float beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

inline complex_double CReflect( Complex<double> alpha )
{ complex_double beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

// Analogues for real variables and integers
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline Int CReflect( Int alpha ) { return alpha; }
inline Int* CReflect( Int* buffer ) { return buffer; }
inline const Int* CReflect( const Int* buffer ) { return buffer; }
/*
inline ElInt CReflect( Int alpha ) { return alpha; }

inline ElInt* CReflect( Int*   buffer ) { return buffer; }
inline Int*   CReflect( ElInt* buffer ) { return buffer; }

inline const ElInt* CReflect( const Int*   buffer ) { return buffer; }
inline const Int*   CReflect( const ElInt* buffer ) { return buffer; }
*/

inline float CReflect( float alpha) { return alpha; }
inline double CReflect( double alpha ) { return alpha; }

inline float* CReflect( float* buffer ) { return buffer; }
inline double* CReflect( double* buffer ) { return buffer; }

inline const float* CReflect( const float* buffer ) { return buffer; }
inline const double* CReflect( const double* buffer ) { return buffer; }

inline ValueInt<Int> CReflect( ElValueInt_i entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_i CReflect( ValueInt<Int> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<Int>* CReflect( ElValueInt_i* entryC )
{ return EL_RC(ValueInt<Int>*,entryC); }
inline ElValueInt_i* CReflect( ValueInt<Int>* entryC )
{ return EL_RC(ElValueInt_i*,entryC); }

inline ValueInt<float> CReflect( ElValueInt_s entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_s CReflect( ValueInt<float> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<float>* CReflect( ElValueInt_s* entryC )
{ return EL_RC(ValueInt<float>*,entryC); }
inline ElValueInt_s* CReflect( ValueInt<float>* entryC )
{ return EL_RC(ElValueInt_s*,entryC); }

inline ValueInt<double> CReflect( ElValueInt_d entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_d CReflect( ValueInt<double> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<double>* CReflect( ElValueInt_d* entryC )
{ return EL_RC(ValueInt<double>*,entryC); }
inline ElValueInt_d* CReflect( ValueInt<double>* entryC )
{ return EL_RC(ElValueInt_d*,entryC); }

inline ValueInt<Complex<float>> CReflect( ElValueInt_c entryC )
{ return {CReflect(entryC.value),entryC.index}; }
inline ElValueInt_c CReflect( ValueInt<Complex<float>> entry )
{ return {CReflect(entry.value),entry.index}; }

inline ValueInt<Complex<float>>* CReflect( ElValueInt_c* entryC )
{ return EL_RC(ValueInt<Complex<float>>*,entryC); }
inline ElValueInt_c* CReflect( ValueInt<Complex<float>>* entryC )
{ return EL_RC(ElValueInt_c*,entryC); }

inline ValueInt<Complex<double>> CReflect( ElValueInt_z entryC )
{ return {CReflect(entryC.value),entryC.index}; }
inline ElValueInt_z CReflect( ValueInt<Complex<double>> entry )
{ return {CReflect(entry.value),entry.index}; }

inline ValueInt<Complex<double>>* CReflect( ElValueInt_z* entryC )
{ return EL_RC(ValueInt<Complex<double>>*,entryC); }
inline ElValueInt_z* CReflect( ValueInt<Complex<double>>* entryC )
{ return EL_RC(ElValueInt_z*,entryC); }

inline Entry<Int> CReflect( ElEntry_i entryC )
{ return { entryC.i, entryC.j, entryC.value }; }
inline ElEntry_i CReflect( Entry<Int> entry )
{ return { entry.i, entry.j, entry.value }; }

inline Entry<float> CReflect( ElEntry_s entryC )
{ return { entryC.i, entryC.j, entryC.value }; }
inline ElEntry_s CReflect( Entry<float> entry )
{ return { entry.i, entry.j, entry.value }; }

inline Entry<double> CReflect( ElEntry_d entryC )
{ return { entryC.i, entryC.j, entryC.value }; }
inline ElEntry_d CReflect( Entry<double> entry )
{ return { entry.i, entry.j, entry.value }; }

inline Entry<Complex<float>> CReflect( ElEntry_c entryC ) 
{ return { entryC.i, entryC.j, CReflect(entryC.value) }; }
inline ElEntry_c CReflect( Entry<Complex<float>> entry )
{ return { entry.i, entry.j, CReflect(entry.value) }; }

inline Entry<Complex<double>> CReflect( ElEntry_z entryC )
{ return { entryC.i, entryC.j, CReflect(entryC.value) }; }
inline ElEntry_z CReflect( Entry<Complex<double>> entry )
{ return { entry.i, entry.j, CReflect(entry.value) }; }

// Matrix
// ------
inline Matrix<Int>* CReflect( ElMatrix_i A )
{ return EL_RC(Matrix<Int>*,A); }

inline Matrix<float>* CReflect( ElMatrix_s A )
{ return EL_RC(Matrix<float>*,A); }

inline Matrix<double>* CReflect( ElMatrix_d A )
{ return EL_RC(Matrix<double>*,A); }

inline Matrix<Complex<float>>* CReflect( ElMatrix_c A )
{ return EL_RC(Matrix<Complex<float>>*,A); }

inline Matrix<Complex<double>>* CReflect( ElMatrix_z A )
{ return EL_RC(Matrix<Complex<double>>*,A); }

inline const Matrix<Int>* CReflect( ElConstMatrix_i A )
{ return EL_RC(const Matrix<Int>*,A); }

inline const Matrix<float>* CReflect( ElConstMatrix_s A )
{ return EL_RC(const Matrix<float>*,A); }

inline const Matrix<double>* CReflect( ElConstMatrix_d A )
{ return EL_RC(const Matrix<double>*,A); }

inline const Matrix<Complex<float>>* CReflect( ElConstMatrix_c A )
{ return EL_RC(const Matrix<Complex<float>>*,A); }

inline const Matrix<Complex<double>>* CReflect( ElConstMatrix_z A )
{ return EL_RC(const Matrix<Complex<double>>*,A); }

inline ElMatrix_i CReflect( Matrix<Int>* A )
{ return (ElMatrix_i)EL_RC(struct ElMatrix_iDummy*,A); }

inline ElMatrix_s CReflect( Matrix<float>* A )
{ return (ElMatrix_s)EL_RC(struct ElMatrix_sDummy*,A); }

inline ElMatrix_d CReflect( Matrix<double>* A )
{ return (ElMatrix_d)EL_RC(struct ElMatrix_dDummy*,A); }

inline ElMatrix_c CReflect( Matrix<Complex<float>>* A )
{ return (ElMatrix_c)EL_RC(struct ElMatrix_cDummy*,A); }

inline ElMatrix_z CReflect( Matrix<Complex<double>>* A )
{ return (ElMatrix_z)EL_RC(struct ElMatrix_zDummy*,A); }

inline ElConstMatrix_i CReflect( const Matrix<Int>* A )
{ return (ElConstMatrix_i)EL_RC(const struct ElMatrix_iDummy*,A); }

inline ElConstMatrix_s CReflect( const Matrix<float>* A )
{ return (ElConstMatrix_s)EL_RC(const struct ElMatrix_sDummy*,A); }

inline ElConstMatrix_d CReflect( const Matrix<double>* A )
{ return (ElConstMatrix_d)EL_RC(const struct ElMatrix_dDummy*,A); }

inline ElConstMatrix_c CReflect( const Matrix<Complex<float>>* A )
{ return (ElConstMatrix_c)EL_RC(const struct ElMatrix_cDummy*,A); }

inline ElConstMatrix_z CReflect( const Matrix<Complex<double>>* A )
{ return (ElConstMatrix_z)EL_RC(const struct ElMatrix_zDummy*,A); }

// ElementalMatrix
// ---------------
inline ElementalMatrix<Int>* 
CReflect( ElDistMatrix_i A )
{ return EL_RC(ElementalMatrix<Int>*,A); }

inline ElementalMatrix<float>* 
CReflect( ElDistMatrix_s A )
{ return EL_RC(ElementalMatrix<float>*,A); }

inline ElementalMatrix<double>* 
CReflect( ElDistMatrix_d A )
{ return EL_RC(ElementalMatrix<double>*,A); }

inline ElementalMatrix<Complex<float>>* 
CReflect( ElDistMatrix_c A )
{ return EL_RC(ElementalMatrix<Complex<float>>*,A); }

inline ElementalMatrix<Complex<double>>* 
CReflect( ElDistMatrix_z A )
{ return EL_RC(ElementalMatrix<Complex<double>>*,A); }

inline const ElementalMatrix<Int>* 
CReflect( ElConstDistMatrix_i A )
{ return EL_RC(const ElementalMatrix<Int>*,A); }

inline const ElementalMatrix<float>* 
CReflect( ElConstDistMatrix_s A )
{ return EL_RC(const ElementalMatrix<float>*,A); }

inline const ElementalMatrix<double>* 
CReflect( ElConstDistMatrix_d A )
{ return EL_RC(const ElementalMatrix<double>*,A); }

inline const ElementalMatrix<Complex<float>>* 
CReflect( ElConstDistMatrix_c A )
{ return EL_RC(const ElementalMatrix<Complex<float>>*,A); }

inline const ElementalMatrix<Complex<double>>* 
CReflect( ElConstDistMatrix_z A )
{ return EL_RC(const ElementalMatrix<Complex<double>>*,A); }

inline ElDistMatrix_i
CReflect( ElementalMatrix<Int>* A )
{ return (ElDistMatrix_i)EL_RC(struct ElDistMatrix_iDummy*,A); }

inline ElDistMatrix_s 
CReflect( ElementalMatrix<float>* A )
{ return (ElDistMatrix_s)EL_RC(struct ElDistMatrix_sDummy*,A); }

inline ElDistMatrix_d 
CReflect( ElementalMatrix<double>* A )
{ return (ElDistMatrix_d)EL_RC(struct ElDistMatrix_dDummy*,A); }

inline ElDistMatrix_c 
CReflect( ElementalMatrix<Complex<float>>* A )
{ return (ElDistMatrix_c)EL_RC(struct ElDistMatrix_cDummy*,A); }

inline ElDistMatrix_z 
CReflect( ElementalMatrix<Complex<double>>* A )
{ return (ElDistMatrix_z)EL_RC(struct ElDistMatrix_zDummy*,A); }

inline ElConstDistMatrix_i
CReflect( const ElementalMatrix<Int>* A )
{ return (ElConstDistMatrix_i)EL_RC(const struct ElDistMatrix_iDummy*,A); }

inline ElConstDistMatrix_s 
CReflect( const ElementalMatrix<float>* A )
{ return (ElConstDistMatrix_s)EL_RC(const struct ElDistMatrix_sDummy*,A); }

inline ElConstDistMatrix_d 
CReflect( const ElementalMatrix<double>* A )
{ return (ElConstDistMatrix_d)EL_RC(const struct ElDistMatrix_dDummy*,A); }

inline ElConstDistMatrix_c 
CReflect( const ElementalMatrix<Complex<float>>* A )
{ return (ElConstDistMatrix_c)EL_RC(const struct ElDistMatrix_cDummy*,A); }

inline ElConstDistMatrix_z 
CReflect( const ElementalMatrix<Complex<double>>* A )
{ return (ElConstDistMatrix_z)EL_RC(const struct ElDistMatrix_zDummy*,A); }

/* Graph
   ----- */
inline Graph* CReflect( ElGraph graph )
{ return EL_RC(Graph*,graph); }

inline const Graph* CReflect( ElConstGraph graph )
{ return EL_RC(const Graph*,graph); }

inline ElGraph CReflect( Graph* graph )
{ return EL_RC(ElGraph,graph); }

inline ElConstGraph CReflect( const Graph* graph )
{ return EL_RC(ElConstGraph,graph); }

/* DistGraph
   --------- */
inline DistGraph* CReflect( ElDistGraph graph )
{ return EL_RC(DistGraph*,graph); }

inline const DistGraph* CReflect( ElConstDistGraph graph )
{ return EL_RC(const DistGraph*,graph); }

inline ElDistGraph CReflect( DistGraph* graph )
{ return EL_RC(ElDistGraph,graph); }

inline ElConstDistGraph CReflect( const DistGraph* graph )
{ return EL_RC(ElConstDistGraph,graph); }

/* SparseMatrix
   ------------ */
inline SparseMatrix<Int>* CReflect( ElSparseMatrix_i A )
{ return EL_RC(SparseMatrix<Int>*,A); }

inline SparseMatrix<float>* CReflect( ElSparseMatrix_s A )
{ return EL_RC(SparseMatrix<float>*,A); }

inline SparseMatrix<double>* CReflect( ElSparseMatrix_d A )
{ return EL_RC(SparseMatrix<double>*,A); }

inline SparseMatrix<Complex<float>>* CReflect( ElSparseMatrix_c A )
{ return EL_RC(SparseMatrix<Complex<float>>*,A); }

inline SparseMatrix<Complex<double>>* CReflect( ElSparseMatrix_z A )
{ return EL_RC(SparseMatrix<Complex<double>>*,A); }

inline const SparseMatrix<Int>* CReflect( ElConstSparseMatrix_i A )
{ return EL_RC(const SparseMatrix<Int>*,A); }

inline const SparseMatrix<float>* CReflect( ElConstSparseMatrix_s A )
{ return EL_RC(const SparseMatrix<float>*,A); }

inline const SparseMatrix<double>* CReflect( ElConstSparseMatrix_d A )
{ return EL_RC(const SparseMatrix<double>*,A); }

inline const SparseMatrix<Complex<float>>* CReflect( ElConstSparseMatrix_c A )
{ return EL_RC(const SparseMatrix<Complex<float>>*,A); }

inline const SparseMatrix<Complex<double>>* CReflect( ElConstSparseMatrix_z A )
{ return EL_RC(const SparseMatrix<Complex<double>>*,A); }

inline ElSparseMatrix_i CReflect( SparseMatrix<Int>* A )
{ return (ElSparseMatrix_i)EL_RC(struct ElSparseMatrix_iDummy*,A); }

inline ElSparseMatrix_s CReflect( SparseMatrix<float>* A )
{ return (ElSparseMatrix_s)EL_RC(struct ElSparseMatrix_sDummy*,A); }

inline ElSparseMatrix_d CReflect( SparseMatrix<double>* A )
{ return (ElSparseMatrix_d)EL_RC(struct ElSparseMatrix_dDummy*,A); }

inline ElSparseMatrix_c CReflect( SparseMatrix<Complex<float>>* A )
{ return (ElSparseMatrix_c)EL_RC(struct ElSparseMatrix_cDummy*,A); }

inline ElSparseMatrix_z CReflect( SparseMatrix<Complex<double>>* A )
{ return (ElSparseMatrix_z)EL_RC(struct ElSparseMatrix_zDummy*,A); }

inline ElConstSparseMatrix_i CReflect( const SparseMatrix<Int>* A )
{ return (ElConstSparseMatrix_i)EL_RC(const struct ElSparseMatrix_iDummy*,A); }

inline ElConstSparseMatrix_s CReflect( const SparseMatrix<float>* A )
{ return (ElConstSparseMatrix_s)EL_RC(const struct ElSparseMatrix_sDummy*,A); }

inline ElConstSparseMatrix_d CReflect( const SparseMatrix<double>* A )
{ return (ElConstSparseMatrix_d)EL_RC(const struct ElSparseMatrix_dDummy*,A); }

inline ElConstSparseMatrix_c CReflect( const SparseMatrix<Complex<float>>* A )
{ return (ElConstSparseMatrix_c)EL_RC(const struct ElSparseMatrix_cDummy*,A); }

inline ElConstSparseMatrix_z CReflect( const SparseMatrix<Complex<double>>* A )
{ return (ElConstSparseMatrix_z)EL_RC(const struct ElSparseMatrix_zDummy*,A); }

/* DistSparseMatrix
   ---------------- */
inline DistSparseMatrix<Int>* CReflect( ElDistSparseMatrix_i A )
{ return EL_RC(DistSparseMatrix<Int>*,A); }

inline DistSparseMatrix<float>* CReflect( ElDistSparseMatrix_s A )
{ return EL_RC(DistSparseMatrix<float>*,A); }

inline DistSparseMatrix<double>* CReflect( ElDistSparseMatrix_d A )
{ return EL_RC(DistSparseMatrix<double>*,A); }

inline DistSparseMatrix<Complex<float>>* CReflect( ElDistSparseMatrix_c A )
{ return EL_RC(DistSparseMatrix<Complex<float>>*,A); }

inline DistSparseMatrix<Complex<double>>* CReflect( ElDistSparseMatrix_z A )
{ return EL_RC(DistSparseMatrix<Complex<double>>*,A); }

inline const DistSparseMatrix<Int>* CReflect( ElConstDistSparseMatrix_i A )
{ return EL_RC(const DistSparseMatrix<Int>*,A); }

inline const DistSparseMatrix<float>* CReflect( ElConstDistSparseMatrix_s A )
{ return EL_RC(const DistSparseMatrix<float>*,A); }

inline const DistSparseMatrix<double>* CReflect( ElConstDistSparseMatrix_d A )
{ return EL_RC(const DistSparseMatrix<double>*,A); }

inline const DistSparseMatrix<Complex<float>>* CReflect
( ElConstDistSparseMatrix_c A )
{ return EL_RC(const DistSparseMatrix<Complex<float>>*,A); }

inline const DistSparseMatrix<Complex<double>>* CReflect
( ElConstDistSparseMatrix_z A )
{ return EL_RC(const DistSparseMatrix<Complex<double>>*,A); }

inline ElDistSparseMatrix_i CReflect( DistSparseMatrix<Int>* A )
{ return (ElDistSparseMatrix_i)EL_RC(struct ElDistSparseMatrix_iDummy*,A); }

inline ElDistSparseMatrix_s CReflect( DistSparseMatrix<float>* A )
{ return (ElDistSparseMatrix_s)EL_RC(struct ElDistSparseMatrix_sDummy*,A); }

inline ElDistSparseMatrix_d CReflect( DistSparseMatrix<double>* A )
{ return (ElDistSparseMatrix_d)EL_RC(struct ElDistSparseMatrix_dDummy*,A); }

inline ElDistSparseMatrix_c CReflect( DistSparseMatrix<Complex<float>>* A )
{ return (ElDistSparseMatrix_c)EL_RC(struct ElDistSparseMatrix_cDummy*,A); }

inline ElDistSparseMatrix_z CReflect( DistSparseMatrix<Complex<double>>* A )
{ return (ElDistSparseMatrix_z)EL_RC(struct ElDistSparseMatrix_zDummy*,A); }

inline ElConstDistSparseMatrix_i CReflect( const DistSparseMatrix<Int>* A )
{ return (ElConstDistSparseMatrix_i)
  EL_RC(const struct ElDistSparseMatrix_iDummy*,A); }

inline ElConstDistSparseMatrix_s CReflect
( const DistSparseMatrix<float>* A )
{ return (ElConstDistSparseMatrix_s)
  EL_RC(const struct ElDistSparseMatrix_sDummy*,A); }

inline ElConstDistSparseMatrix_d CReflect
( const DistSparseMatrix<double>* A )
{ return (ElConstDistSparseMatrix_d)
  EL_RC(const struct ElDistSparseMatrix_dDummy*,A); }

inline ElConstDistSparseMatrix_c CReflect
( const DistSparseMatrix<Complex<float>>* A )
{ return (ElConstDistSparseMatrix_c)
  EL_RC(const struct ElDistSparseMatrix_cDummy*,A); }

inline ElConstDistSparseMatrix_z CReflect
( const DistSparseMatrix<Complex<double>>* A )
{ return (ElConstDistSparseMatrix_z)
  EL_RC(const struct ElDistSparseMatrix_zDummy*,A); }

/* DistMultiVec
   ------------ */
inline DistMultiVec<Int>* CReflect( ElDistMultiVec_i A )
{ return EL_RC(DistMultiVec<Int>*,A); }

inline DistMultiVec<float>* CReflect( ElDistMultiVec_s A )
{ return EL_RC(DistMultiVec<float>*,A); }

inline DistMultiVec<double>* CReflect( ElDistMultiVec_d A )
{ return EL_RC(DistMultiVec<double>*,A); }

inline DistMultiVec<Complex<float>>* CReflect( ElDistMultiVec_c A )
{ return EL_RC(DistMultiVec<Complex<float>>*,A); }

inline DistMultiVec<Complex<double>>* CReflect( ElDistMultiVec_z A )
{ return EL_RC(DistMultiVec<Complex<double>>*,A); }

inline const DistMultiVec<Int>* CReflect( ElConstDistMultiVec_i A )
{ return EL_RC(const DistMultiVec<Int>*,A); }

inline const DistMultiVec<float>* CReflect( ElConstDistMultiVec_s A )
{ return EL_RC(const DistMultiVec<float>*,A); }

inline const DistMultiVec<double>* CReflect( ElConstDistMultiVec_d A )
{ return EL_RC(const DistMultiVec<double>*,A); }

inline const DistMultiVec<Complex<float>>* CReflect( ElConstDistMultiVec_c A )
{ return EL_RC(const DistMultiVec<Complex<float>>*,A); }

inline const DistMultiVec<Complex<double>>* CReflect( ElConstDistMultiVec_z A )
{ return EL_RC(const DistMultiVec<Complex<double>>*,A); }

inline ElDistMultiVec_i CReflect( DistMultiVec<Int>* A )
{ return (ElDistMultiVec_i)EL_RC(struct ElDistMultiVec_iDummy*,A); }

inline ElDistMultiVec_s CReflect( DistMultiVec<float>* A )
{ return (ElDistMultiVec_s)EL_RC(struct ElDistMultiVec_sDummy*,A); }

inline ElDistMultiVec_d CReflect( DistMultiVec<double>* A )
{ return (ElDistMultiVec_d)EL_RC(struct ElDistMultiVec_dDummy*,A); }

inline ElDistMultiVec_c CReflect( DistMultiVec<Complex<float>>* A )
{ return (ElDistMultiVec_c)EL_RC(struct ElDistMultiVec_cDummy*,A); }

inline ElDistMultiVec_z CReflect( DistMultiVec<Complex<double>>* A )
{ return (ElDistMultiVec_z)EL_RC(struct ElDistMultiVec_zDummy*,A); }

inline ElConstDistMultiVec_i CReflect( const DistMultiVec<Int>* A )
{ return (ElConstDistMultiVec_i)EL_RC(const struct ElDistMultiVec_iDummy*,A); }

inline ElConstDistMultiVec_s CReflect( const DistMultiVec<float>* A )
{ return (ElConstDistMultiVec_s)EL_RC(const struct ElDistMultiVec_sDummy*,A); }

inline ElConstDistMultiVec_d CReflect( const DistMultiVec<double>* A )
{ return (ElConstDistMultiVec_d)EL_RC(const struct ElDistMultiVec_dDummy*,A); }

inline ElConstDistMultiVec_c CReflect( const DistMultiVec<Complex<float>>* A )
{ return (ElConstDistMultiVec_c)EL_RC(const struct ElDistMultiVec_cDummy*,A); }

inline ElConstDistMultiVec_z CReflect( const DistMultiVec<Complex<double>>* A )
{ return (ElConstDistMultiVec_z)EL_RC(const struct ElDistMultiVec_zDummy*,A); }

inline ElDistData CReflect( const DistData& data )
{
    ElDistData distData;
    distData.colDist = CReflect(data.colDist);
    distData.rowDist = CReflect(data.rowDist);
    distData.blockHeight = data.blockHeight;
    distData.blockWidth = data.blockWidth;
    distData.colAlign = data.colAlign;
    distData.rowAlign = data.rowAlign;
    distData.colCut = data.colCut;
    distData.rowCut = data.rowCut;
    distData.root = data.root;
    distData.grid = CReflect(data.grid);
    return distData;
}

inline DistData CReflect( const ElDistData& distData )
{
    DistData data;
    data.colDist = CReflect(distData.colDist);
    data.rowDist = CReflect(distData.rowDist);
    data.blockHeight = distData.blockHeight;
    data.blockWidth = distData.blockWidth;
    data.colAlign = distData.colAlign;
    data.rowAlign = distData.rowAlign;
    data.colCut = distData.colCut;
    data.rowCut = distData.rowCut;
    data.root = distData.root;
    data.grid = CReflect(distData.grid);
    return data;
}

inline ElSafeProduct_s CReflect( const SafeProduct<float>& prod )
{ 
    ElSafeProduct_s prodC;    
    prodC.rho = prod.rho;
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_d CReflect( const SafeProduct<double>& prod )
{ 
    ElSafeProduct_d prodC;    
    prodC.rho = prod.rho;
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_c CReflect( const SafeProduct<Complex<float>>& prod )
{ 
    ElSafeProduct_c prodC;    
    prodC.rho = CReflect(prod.rho);
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}
inline ElSafeProduct_z CReflect( const SafeProduct<Complex<double>>& prod )
{ 
    ElSafeProduct_z prodC;    
    prodC.rho = CReflect(prod.rho);
    prodC.kappa = prod.kappa;
    prodC.n = prod.n;
    return prodC;
}

inline SafeProduct<float> CReflect( const ElSafeProduct_s& prodC )
{ 
    SafeProduct<float> prod( prodC.n );
    prod.rho = prodC.rho;
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<double> CReflect( const ElSafeProduct_d& prodC )
{ 
    SafeProduct<double> prod( prodC.n );
    prod.rho = prodC.rho;
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<Complex<float>> CReflect( const ElSafeProduct_c& prodC )
{ 
    SafeProduct<Complex<float>> prod( prodC.n );
    prod.rho = CReflect(prodC.rho);
    prod.kappa = prodC.kappa;
    return prod;
}
inline SafeProduct<Complex<double>> CReflect( const ElSafeProduct_z& prodC )
{ 
    SafeProduct<Complex<double>> prod( prodC.n );
    prod.rho = CReflect(prodC.rho);
    prod.kappa = prodC.kappa;
    return prod;
}

inline ElFileFormat CReflect( FileFormat format )
{ return static_cast<ElFileFormat>(format); }
inline FileFormat CReflect( ElFileFormat format )
{ return static_cast<FileFormat>(format); }

inline ElPermutationMeta CReflect( const PermutationMeta& meta )
{
    ElPermutationMeta metaC;

    metaC.align = meta.align;
    metaC.comm = meta.comm.comm;

    const Int commSize = mpi::Size( meta.comm );
    metaC.sendCounts = new int[commSize];
    metaC.sendDispls = new int[commSize];
    metaC.recvCounts = new int[commSize];
    metaC.recvDispls = new int[commSize];
    MemCopy( metaC.sendCounts, meta.sendCounts.data(), commSize );
    MemCopy( metaC.sendDispls, meta.sendDispls.data(), commSize );
    MemCopy( metaC.recvCounts, meta.recvCounts.data(), commSize );
    MemCopy( metaC.recvDispls, meta.recvDispls.data(), commSize );

    metaC.numSendIdx = meta.sendIdx.size();
    metaC.numRecvIdx = meta.recvIdx.size();
    metaC.sendIdx   = new int[metaC.numSendIdx];
    metaC.sendRanks = new int[metaC.numSendIdx];
    metaC.recvIdx   = new int[metaC.numRecvIdx];
    metaC.recvRanks = new int[metaC.numRecvIdx];
    MemCopy( metaC.sendIdx,   meta.sendIdx.data(),   metaC.numSendIdx );
    MemCopy( metaC.sendRanks, meta.sendRanks.data(), metaC.numSendIdx );
    MemCopy( metaC.recvIdx,   meta.recvIdx.data(),   metaC.numRecvIdx );
    MemCopy( metaC.recvRanks, meta.recvRanks.data(), metaC.numRecvIdx );

    return metaC;
}

inline PermutationMeta CReflect( const ElPermutationMeta& metaC )
{
    PermutationMeta meta;

    meta.align = metaC.align;
    meta.comm = metaC.comm;

    int commSize;
    MPI_Comm_size( metaC.comm, &commSize );
    meta.sendCounts =
        vector<int>( metaC.sendCounts, metaC.sendCounts+commSize );
    meta.sendDispls =
        vector<int>( metaC.sendDispls, metaC.sendDispls+commSize );
    meta.recvCounts =
        vector<int>( metaC.recvCounts, metaC.recvCounts+commSize );
    meta.recvDispls =
        vector<int>( metaC.recvDispls, metaC.recvDispls+commSize );

    meta.sendIdx =
        vector<int>( metaC.sendIdx, metaC.sendIdx+metaC.numSendIdx );
    meta.sendRanks =
        vector<int>( metaC.sendRanks, metaC.sendRanks+metaC.numSendIdx );
    meta.recvIdx =
        vector<int>( metaC.recvIdx, metaC.recvIdx+metaC.numRecvIdx );
    meta.recvRanks =
        vector<int>( metaC.recvRanks, metaC.recvRanks+metaC.numRecvIdx );

    return meta;
}

inline Permutation* CReflect( ElPermutation p )
{ return EL_RC(Permutation*,p); }
inline ElPermutation CReflect( Permutation* p )
{ return (ElPermutation)EL_RC(struct ElPermutationDummy*,p); }

inline DistPermutation* CReflect( ElDistPermutation p )
{ return EL_RC(DistPermutation*,p); }
inline ElDistPermutation CReflect( DistPermutation* p )
{ return (ElDistPermutation)EL_RC(struct ElDistPermutationDummy*,p); }

inline const Permutation* CReflect( ElConstPermutation p )
{ return EL_RC(const Permutation*,p); }
inline ElConstPermutation CReflect( const Permutation* p )
{ return (ElConstPermutation)EL_RC(const struct ElPermutationDummy*,p); }

inline const DistPermutation* CReflect( ElConstDistPermutation p )
{ return EL_RC(const DistPermutation*,p); }
inline ElConstDistPermutation CReflect( const DistPermutation* p )
{ return (ElConstDistPermutation)
         EL_RC(const struct ElDistPermutationDummy*,p); }

} // namespace El

#endif // ifndef EL_CORE_CREFLECT_C_HPP
