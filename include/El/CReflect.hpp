/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CREFLECT_C_HPP
#define EL_CREFLECT_C_HPP

// TODO: Break this file into pieces. It is getting ridiculous.

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

inline ElElementalData CReflect( const ElementalData& data )
{
    ElElementalData distData;
    distData.colDist = CReflect(data.colDist);
    distData.rowDist = CReflect(data.rowDist);
    distData.colAlign = data.colAlign;
    distData.rowAlign = data.rowAlign;
    distData.root = data.root;
    distData.grid = CReflect(data.grid);
    return distData;
}

inline ElementalData CReflect( const ElElementalData& distData )
{
    ElementalData data;
    data.colDist = CReflect(distData.colDist);
    data.rowDist = CReflect(distData.rowDist);
    data.colAlign = distData.colAlign;
    data.rowAlign = distData.rowAlign;
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

// Input/Output
// ------------
inline ElFileFormat CReflect( FileFormat format )
{ return static_cast<ElFileFormat>(format); }
inline FileFormat CReflect( ElFileFormat format )
{ return static_cast<FileFormat>(format); }

inline ElColorMap CReflect( ColorMap map )
{ return static_cast<ElColorMap>(map); }
inline ColorMap CReflect( ElColorMap map )
{ return static_cast<ColorMap>(map); }

// BLAS-like
// ---------
inline ElGemmAlgorithm CReflect( GemmAlgorithm alg )
{ return static_cast<ElGemmAlgorithm>(alg); }
inline GemmAlgorithm CReflect( ElGemmAlgorithm alg )
{ return static_cast<GemmAlgorithm>(alg); }

template<typename T>
inline ElSymvCtrl
CReflect( const SymvCtrl<T>& ctrl )
{ 
    ElSymvCtrl ctrlC;
    ctrlC.bsize = ctrl.bsize;
    ctrlC.avoidTrmvBasedLocalSymv = ctrl.avoidTrmvBasedLocalSymv;
    return ctrlC;
}

template<typename T>
inline SymvCtrl<T>
CReflect( const ElSymvCtrl& ctrlC )
{ 
    SymvCtrl<T> ctrl;
    ctrl.bsize = ctrlC.bsize;
    ctrl.avoidTrmvBasedLocalSymv = ctrlC.avoidTrmvBasedLocalSymv;
    return ctrl;
}

// LAPACK-like
// -----------

inline ElSortType CReflect( SortType type )
{ return static_cast<ElSortType>(type); }

inline SortType CReflect( ElSortType type )
{ return static_cast<SortType>(type); }

// Permutations
// ^^^^^^^^^^^^

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

// Condensed form
// ^^^^^^^^^^^^^^
inline ElHermitianTridiagApproach 
CReflect( HermitianTridiagApproach approach )
{ return static_cast<ElHermitianTridiagApproach>( approach ); }

inline HermitianTridiagApproach 
CReflect( ElHermitianTridiagApproach approach )
{ return static_cast<HermitianTridiagApproach>( approach ); }

template<typename F>
inline ElHermitianTridiagCtrl
CReflect( const HermitianTridiagCtrl<F>& ctrl )
{ 
    ElHermitianTridiagCtrl ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.order = CReflect(ctrl.order);
    ctrlC.symvCtrl = CReflect(ctrl.symvCtrl);
    return ctrlC;
}

template<typename F>
inline HermitianTridiagCtrl<F>
CReflect( const ElHermitianTridiagCtrl& ctrlC )
{ 
    HermitianTridiagCtrl<F> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.order = CReflect(ctrlC.order);
    ctrl.symvCtrl = CReflect<F>(ctrlC.symvCtrl);
    return ctrl;
}

// Decompositions
// ^^^^^^^^^^^^^^

/* Pencil */
inline ElPencil CReflect( Pencil pencil )
{ return static_cast<ElPencil>(pencil); }

inline Pencil CReflect( ElPencil pencil )
{ return static_cast<Pencil>(pencil); }

/* HermitianSDCCtrl */
inline ElHermitianSDCCtrl_s CReflect( const HermitianSDCCtrl<float>& ctrl )
{
    ElHermitianSDCCtrl_s ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElHermitianSDCCtrl_d CReflect( const HermitianSDCCtrl<double>& ctrl )
{
    ElHermitianSDCCtrl_d ctrlC;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline HermitianSDCCtrl<float> CReflect( const ElHermitianSDCCtrl_s& ctrlC )
{
    HermitianSDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline HermitianSDCCtrl<double> CReflect( const ElHermitianSDCCtrl_d& ctrlC )
{
    HermitianSDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* HermitianEigSubset */
inline ElHermitianEigSubset_s CReflect
( const HermitianEigSubset<float>& subset )
{
    ElHermitianEigSubset_s subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}
inline ElHermitianEigSubset_d CReflect
( const HermitianEigSubset<double>& subset )
{
    ElHermitianEigSubset_d subsetC;
    subsetC.indexSubset = subset.indexSubset;
    subsetC.lowerIndex = subset.lowerIndex;
    subsetC.upperIndex = subset.upperIndex;
    subsetC.rangeSubset = subset.rangeSubset;
    subsetC.lowerBound = subset.lowerBound;
    subsetC.upperBound = subset.upperBound;
    return subsetC;
}

inline HermitianEigSubset<float> CReflect
( const ElHermitianEigSubset_s& subsetC )
{
    HermitianEigSubset<float> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}
inline HermitianEigSubset<double> CReflect
( const ElHermitianEigSubset_d& subsetC )
{
    HermitianEigSubset<double> subset;
    subset.indexSubset = subsetC.indexSubset;
    subset.lowerIndex = subsetC.lowerIndex;
    subset.upperIndex = subsetC.upperIndex;
    subset.rangeSubset = subsetC.rangeSubset;
    subset.lowerBound = subsetC.lowerBound;
    subset.upperBound = subsetC.upperBound;
    return subset;
}

/* HermitianEigCtrl */
inline ElHermitianEigCtrl_s CReflect( const HermitianEigCtrl<float>& ctrl )
{
    ElHermitianEigCtrl_s ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline ElHermitianEigCtrl_d CReflect( const HermitianEigCtrl<double>& ctrl )
{
    ElHermitianEigCtrl_d ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_c 
CReflect( const HermitianEigCtrl<Complex<float>>& ctrl )
{
    ElHermitianEigCtrl_c ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}
inline ElHermitianEigCtrl_z
CReflect( const HermitianEigCtrl<Complex<double>>& ctrl )
{
    ElHermitianEigCtrl_z ctrlC;
    ctrlC.tridiagCtrl = CReflect( ctrl.tridiagCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    ctrlC.useSDC = ctrl.useSDC;
    return ctrlC;
}

inline HermitianEigCtrl<float> CReflect( const ElHermitianEigCtrl_s& ctrlC )
{
    HermitianEigCtrl<float> ctrl;
    ctrl.tridiagCtrl = CReflect<float>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<double> CReflect( const ElHermitianEigCtrl_d& ctrlC )
{
    HermitianEigCtrl<double> ctrl;
    ctrl.tridiagCtrl = CReflect<double>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<float>> 
CReflect( const ElHermitianEigCtrl_c& ctrlC )
{
    HermitianEigCtrl<Complex<float>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<float>>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}
inline HermitianEigCtrl<Complex<double>> 
CReflect( const ElHermitianEigCtrl_z& ctrlC )
{
    HermitianEigCtrl<Complex<double>> ctrl;
    ctrl.tridiagCtrl = CReflect<Complex<double>>( ctrlC.tridiagCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    ctrl.useSDC = ctrlC.useSDC;
    return ctrl;
}

/* QDWHCtrl */
inline ElQDWHCtrl CReflect( const QDWHCtrl& ctrl )
{
    ElQDWHCtrl ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.maxIts = ctrl.maxIts;
    return ctrlC;
}

inline QDWHCtrl CReflect( const ElQDWHCtrl& ctrlC )
{
    QDWHCtrl ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.maxIts = ctrlC.maxIts;
    return ctrl;
}

/* PolarCtrl */
inline ElPolarCtrl CReflect( const PolarCtrl& ctrl )
{
    ElPolarCtrl ctrlC;
    ctrlC.qdwh = ctrl.qdwh;
    ctrlC.qdwhCtrl = CReflect(ctrl.qdwhCtrl);
    return ctrlC;
}

inline PolarCtrl CReflect( const ElPolarCtrl& ctrlC )
{
    PolarCtrl ctrl;
    ctrl.qdwh = ctrlC.qdwh;
    ctrl.qdwhCtrl = CReflect(ctrlC.qdwhCtrl);
    return ctrl;
}

/* QDWHInfo */
inline ElQDWHInfo CReflect( const QDWHInfo& info )
{
    ElQDWHInfo infoC;
    infoC.numIts = info.numIts;
    infoC.numQRIts = info.numQRIts;
    infoC.numCholIts = info.numCholIts;
    return infoC;
}

inline QDWHInfo CReflect( const ElQDWHInfo& infoC )
{
    QDWHInfo info;
    info.numIts = infoC.numIts;
    info.numQRIts = infoC.numQRIts;
    info.numCholIts = infoC.numCholIts;
    return info;
}

/* PolarInfo */
inline ElPolarInfo CReflect( const PolarInfo& info )
{
    ElPolarInfo infoC;
    infoC.qdwhInfo = CReflect(info.qdwhInfo);
    return infoC;
}

inline PolarInfo CReflect( const ElPolarInfo& infoC )
{
    PolarInfo info;
    info.qdwhInfo = CReflect(infoC.qdwhInfo);
    return info;
}

/* SVDCtrl */
inline ElSVDApproach CReflect( SVDApproach approach )
{ return static_cast<ElSVDApproach>( approach ); }

inline SVDApproach CReflect( ElSVDApproach approach )
{ return static_cast<SVDApproach>( approach ); }

inline SVDCtrl<float> CReflect( const ElSVDCtrl_s& ctrlC )
{
    SVDCtrl<float> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.avoidComputingU = ctrlC.avoidComputingU;
    ctrl.avoidComputingV = ctrlC.avoidComputingV;
    ctrl.time = ctrlC.time;
    ctrl.avoidLibflame = ctrlC.avoidLibflame;

    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;

    return ctrl;
}

inline SVDCtrl<double> CReflect( const ElSVDCtrl_d& ctrlC )
{
    SVDCtrl<double> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.overwrite = ctrlC.overwrite;
    ctrl.avoidComputingU = ctrlC.avoidComputingU;
    ctrl.avoidComputingV = ctrlC.avoidComputingV;
    ctrl.time = ctrlC.time;
    ctrl.avoidLibflame = ctrlC.avoidLibflame;

    ctrl.seqQR = ctrlC.seqQR;
    ctrl.valChanRatio = ctrlC.valChanRatio;
    ctrl.fullChanRatio = ctrlC.fullChanRatio;
    ctrl.relative = ctrlC.relative;
    ctrl.tol = ctrlC.tol;

    return ctrl;
}

inline ElSVDCtrl_s CReflect( const SVDCtrl<float>& ctrl )
{
    ElSVDCtrl_s ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.avoidComputingU = ctrl.avoidComputingU;
    ctrlC.avoidComputingV = ctrl.avoidComputingV;
    ctrlC.time = ctrl.time;
    ctrlC.avoidLibflame = ctrl.avoidLibflame;

    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;

    return ctrlC;
}

inline ElSVDCtrl_d CReflect( const SVDCtrl<double>& ctrl )
{
    ElSVDCtrl_d ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.overwrite = ctrl.overwrite;
    ctrlC.avoidComputingU = ctrl.avoidComputingU;
    ctrlC.avoidComputingV = ctrl.avoidComputingV;
    ctrlC.time = ctrl.time;
    ctrlC.avoidLibflame = ctrl.avoidLibflame;

    ctrlC.seqQR = ctrl.seqQR;
    ctrlC.valChanRatio = ctrl.valChanRatio;
    ctrlC.fullChanRatio = ctrl.fullChanRatio;
    ctrlC.relative = ctrl.relative;
    ctrlC.tol = ctrl.tol;

    return ctrlC;
}

/* HessQRCtrl */
inline ElHessQRCtrl CReflect( const HessQRCtrl& ctrl )
{
    ElHessQRCtrl ctrlC;
    ctrlC.distAED = ctrl.distAED;
    ctrlC.blockHeight = ctrl.blockHeight;
    ctrlC.blockWidth = ctrl.blockWidth;
    return ctrlC;
}

inline HessQRCtrl CReflect( const ElHessQRCtrl& ctrlC )
{
    HessQRCtrl ctrl;
    ctrl.distAED = ctrlC.distAED;
    ctrl.blockHeight = ctrlC.blockHeight;
    ctrl.blockWidth = ctrlC.blockWidth;
    return ctrl;
}

/* SDCCtrl */
inline ElSDCCtrl_s CReflect( const SDCCtrl<float>& ctrl )
{
    ElSDCCtrl_s ctrlC;    
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}
inline ElSDCCtrl_d CReflect( const SDCCtrl<double>& ctrl )
{
    ElSDCCtrl_d ctrlC;    
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.maxInnerIts = ctrl.maxInnerIts;
    ctrlC.maxOuterIts = ctrl.maxOuterIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.spreadFactor = ctrl.spreadFactor;
    ctrlC.random = ctrl.random;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline SDCCtrl<float> CReflect( const ElSDCCtrl_s& ctrlC )
{
    SDCCtrl<float> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}
inline SDCCtrl<double> CReflect( const ElSDCCtrl_d& ctrlC )
{
    SDCCtrl<double> ctrl;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.maxInnerIts = ctrlC.maxInnerIts;
    ctrl.maxOuterIts = ctrlC.maxOuterIts;
    ctrl.tol = ctrlC.tol;
    ctrl.spreadFactor = ctrlC.spreadFactor;
    ctrl.random = ctrlC.random;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* SchurCtrl */
inline ElSchurCtrl_s CReflect( const SchurCtrl<float>& ctrl )
{
    ElSchurCtrl_s ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    return ctrlC;
}
inline ElSchurCtrl_d CReflect( const SchurCtrl<double>& ctrl )
{
    ElSchurCtrl_d ctrlC;
    ctrlC.useSDC = ctrl.useSDC;
    ctrlC.qrCtrl = CReflect( ctrl.qrCtrl );
    ctrlC.sdcCtrl = CReflect( ctrl.sdcCtrl );
    return ctrlC;
}

inline SchurCtrl<float> CReflect( const ElSchurCtrl_s& ctrlC )
{
    SchurCtrl<float> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    return ctrl;
}
inline SchurCtrl<double> CReflect( const ElSchurCtrl_d& ctrlC )
{
    SchurCtrl<double> ctrl;
    ctrl.useSDC = ctrlC.useSDC;
    ctrl.qrCtrl = CReflect( ctrlC.qrCtrl );
    ctrl.sdcCtrl = CReflect( ctrlC.sdcCtrl );
    return ctrl;
}

// Factorizations
// ^^^^^^^^^^^^^^
inline ElLDLPivotType CReflect( LDLPivotType pivotType )
{ return static_cast<ElLDLPivotType>( pivotType ); }

inline LDLPivotType CReflect( ElLDLPivotType pivotType )
{ return static_cast<LDLPivotType>( pivotType ); }

inline ElLDLPivotCtrl_s CReflect( const LDLPivotCtrl<float>& ctrl )
{
    ElLDLPivotCtrl_s ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}
inline ElLDLPivotCtrl_d CReflect( const LDLPivotCtrl<double>& ctrl )
{
    ElLDLPivotCtrl_d ctrlC;
    ctrlC.pivotType = CReflect(ctrl.pivotType);
    ctrlC.gamma = ctrl.gamma;
    return ctrlC;
}

inline LDLPivotCtrl<float> CReflect( const ElLDLPivotCtrl_s& ctrlC )
{
    LDLPivotCtrl<float> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}
inline LDLPivotCtrl<double> CReflect( const ElLDLPivotCtrl_d& ctrlC )
{
    LDLPivotCtrl<double> ctrl;
    ctrl.pivotType = CReflect(ctrlC.pivotType);
    ctrl.gamma = ctrlC.gamma;
    return ctrl;
}

inline ElLDLPivot CReflect( const LDLPivot& pivot )
{
    ElLDLPivot pivotC;
    pivotC.nb = pivot.nb;
    pivotC.from[0] = pivot.from[0];
    pivotC.from[1] = pivot.from[1];
    return pivotC;
}

inline LDLPivot CReflect( const ElLDLPivot& pivotC )
{
    LDLPivot pivot;
    pivot.nb = pivotC.nb;
    pivot.from[0] = pivotC.from[0];
    pivot.from[1] = pivotC.from[1];
    return pivot;
}

inline ElInertiaType CReflect( const InertiaType& inertia )
{ 
    ElInertiaType inertiaC;
    inertiaC.numPositive = inertia.numPositive;
    inertiaC.numNegative = inertia.numNegative;
    inertiaC.numZero = inertia.numZero;
    return inertiaC;
}

inline InertiaType CReflect( const ElInertiaType& inertiaC )
{ 
    InertiaType inertia;
    inertia.numPositive = inertiaC.numPositive;
    inertia.numNegative = inertiaC.numNegative;
    inertia.numZero = inertiaC.numZero;
    return inertia;
}

inline ElQRCtrl_s CReflect( const QRCtrl<float>& ctrl )
{ 
    ElQRCtrl_s ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    return ctrlC;
}
inline ElQRCtrl_d CReflect( const QRCtrl<double>& ctrl )
{ 
    ElQRCtrl_d ctrlC;
    ctrlC.colPiv = ctrl.colPiv;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    return ctrlC;
}

inline QRCtrl<float> CReflect( const ElQRCtrl_s& ctrlC )
{ 
    QRCtrl<float> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    return ctrl;
}
inline QRCtrl<double> CReflect( const ElQRCtrl_d& ctrlC )
{ 
    QRCtrl<double> ctrl;
    ctrl.colPiv = ctrlC.colPiv;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    return ctrl;
}

inline RegSolveAlg CReflect( ElRegSolveAlg alg ) 
{ return static_cast<RegSolveAlg>(alg); }
inline ElRegSolveAlg CReflect( RegSolveAlg alg )
{ return static_cast<ElRegSolveAlg>(alg); }

inline ElRegSolveCtrl_s CReflect( const RegSolveCtrl<float>& ctrl )
{
    ElRegSolveCtrl_s ctrlC;
    ctrlC.alg          = CReflect(ctrl.alg);
    ctrlC.relTol       = ctrl.relTol;
    ctrlC.relTolRefine = ctrl.relTolRefine;
    ctrlC.maxRefineIts = ctrl.maxRefineIts;
    ctrlC.restart      = ctrl.restart;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline ElRegSolveCtrl_d CReflect( const RegSolveCtrl<double>& ctrl )
{
    ElRegSolveCtrl_d ctrlC;
    ctrlC.alg          = CReflect(ctrl.alg);
    ctrlC.relTol       = ctrl.relTol;
    ctrlC.relTolRefine = ctrl.relTolRefine;
    ctrlC.maxRefineIts = ctrl.maxRefineIts;
    ctrlC.restart      = ctrl.restart;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline RegSolveCtrl<float> CReflect( const ElRegSolveCtrl_s& ctrlC )
{
    RegSolveCtrl<float> ctrl;
    ctrl.alg          = CReflect(ctrlC.alg);
    ctrl.relTol       = ctrlC.relTol;
    ctrl.relTolRefine = ctrlC.relTolRefine;
    ctrl.maxRefineIts = ctrlC.maxRefineIts;
    ctrl.restart      = ctrlC.restart;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline RegSolveCtrl<double> CReflect( const ElRegSolveCtrl_d& ctrlC )
{
    RegSolveCtrl<double> ctrl;
    ctrl.alg          = CReflect(ctrlC.alg);
    ctrl.relTol       = ctrlC.relTol;
    ctrl.relTolRefine = ctrlC.relTolRefine;
    ctrl.maxRefineIts = ctrlC.maxRefineIts;
    ctrl.restart      = ctrlC.restart;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline ElLeastSquaresCtrl_s CReflect( const LeastSquaresCtrl<float>& ctrl )
{
    ElLeastSquaresCtrl_s ctrlC;
    ctrlC.scaleTwoNorm = ctrl.scaleTwoNorm;
    ctrlC.basisSize    = ctrl.basisSize;
    ctrlC.alpha        = ctrl.alpha; 
    ctrlC.solveCtrl    = CReflect(ctrl.solveCtrl);
    ctrlC.equilibrate  = ctrl.equilibrate;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline ElLeastSquaresCtrl_d CReflect( const LeastSquaresCtrl<double>& ctrl )
{
    ElLeastSquaresCtrl_d ctrlC;
    ctrlC.scaleTwoNorm = ctrl.scaleTwoNorm;
    ctrlC.basisSize    = ctrl.basisSize;
    ctrlC.alpha        = ctrl.alpha; 
    ctrlC.solveCtrl    = CReflect(ctrl.solveCtrl);
    ctrlC.equilibrate  = ctrl.equilibrate;
    ctrlC.progress     = ctrl.progress;
    ctrlC.time         = ctrl.time;
    return ctrlC;
}

inline LeastSquaresCtrl<float> CReflect( const ElLeastSquaresCtrl_s& ctrlC )
{
    LeastSquaresCtrl<float> ctrl;
    ctrl.scaleTwoNorm = ctrlC.scaleTwoNorm;
    ctrl.basisSize    = ctrlC.basisSize;
    ctrl.alpha        = ctrlC.alpha; 
    ctrl.solveCtrl    = CReflect(ctrlC.solveCtrl);
    ctrl.equilibrate  = ctrlC.equilibrate;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

inline LeastSquaresCtrl<double> CReflect( const ElLeastSquaresCtrl_d& ctrlC )
{
    LeastSquaresCtrl<double> ctrl;
    ctrl.scaleTwoNorm = ctrlC.scaleTwoNorm;
    ctrl.basisSize    = ctrlC.basisSize;
    ctrl.alpha        = ctrlC.alpha; 
    ctrl.solveCtrl    = CReflect(ctrlC.solveCtrl);
    ctrl.equilibrate  = ctrlC.equilibrate;
    ctrl.progress     = ctrlC.progress;
    ctrl.time         = ctrlC.time;
    return ctrl;
}

// Properties
// ^^^^^^^^^^
inline ElNormType CReflect( NormType type )
{ return static_cast<ElNormType>(type); }
inline NormType CReflect( ElNormType type )
{ return static_cast<NormType>(type); }

inline ElPseudospecNorm CReflect( PseudospecNorm psNorm )
{ return static_cast<ElPseudospecNorm>(psNorm); }
inline PseudospecNorm CReflect( ElPseudospecNorm psNorm )
{ return static_cast<PseudospecNorm>(psNorm); }

inline ElSnapshotCtrl CReflect( const SnapshotCtrl& ctrl )
{
    ElSnapshotCtrl ctrlC;
    ctrlC.realSize = ctrl.realSize; 
    ctrlC.imagSize = ctrl.imagSize;
    ctrlC.imgSaveFreq = ctrl.imgSaveFreq;
    ctrlC.numSaveFreq = ctrl.numSaveFreq;
    ctrlC.imgDispFreq = ctrl.imgDispFreq;
    ctrlC.imgSaveCount = ctrl.imgSaveCount;
    ctrlC.numSaveCount = ctrl.numSaveCount;
    ctrlC.imgDispCount = ctrl.imgDispCount;
    ctrlC.imgBase = CReflect(ctrl.imgBase);
    ctrlC.numBase = CReflect(ctrl.numBase);
    ctrlC.imgFormat = CReflect(ctrl.imgFormat);
    ctrlC.numFormat = CReflect(ctrl.numFormat);
    ctrlC.itCounts = ctrl.itCounts;
    return ctrlC;
}
inline SnapshotCtrl CReflect( const ElSnapshotCtrl& ctrlC )
{
    SnapshotCtrl ctrl;
    ctrl.realSize = ctrlC.realSize; 
    ctrl.imagSize = ctrlC.imagSize;
    ctrl.imgSaveFreq = ctrlC.imgSaveFreq;
    ctrl.numSaveFreq = ctrlC.numSaveFreq;
    ctrl.imgDispFreq = ctrlC.imgDispFreq;
    ctrl.imgSaveCount = ctrlC.imgSaveCount;
    ctrl.numSaveCount = ctrlC.numSaveCount;
    ctrl.imgDispCount = ctrlC.imgDispCount;
    ctrl.imgBase = CReflect(ctrlC.imgBase);
    ctrl.numBase = CReflect(ctrlC.numBase);
    ctrl.imgFormat = CReflect(ctrlC.imgFormat);
    ctrl.numFormat = CReflect(ctrlC.numFormat);
    ctrl.itCounts = ctrlC.itCounts;
    return ctrl;
}

inline ElPseudospecCtrl_s CReflect( const PseudospecCtrl<float>& ctrl )
{
    ElPseudospecCtrl_s ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}
inline ElPseudospecCtrl_d CReflect( const PseudospecCtrl<double>& ctrl )
{
    ElPseudospecCtrl_d ctrlC;
    ctrlC.norm = CReflect(ctrl.norm);
    ctrlC.blockWidth = ctrl.norm;
    ctrlC.schur = ctrl.schur;
    ctrlC.forceComplexSchur = ctrl.forceComplexSchur;
    ctrlC.forceComplexPs = ctrl.forceComplexPs;
    ctrlC.schurCtrl = CReflect(ctrl.schurCtrl);
    ctrlC.maxIts = ctrl.maxIts;
    ctrlC.tol = ctrl.tol;
    ctrlC.deflate = ctrl.deflate;
    ctrlC.arnoldi = ctrl.arnoldi;
    ctrlC.basisSize = ctrl.basisSize;
    ctrlC.reorthog = ctrl.reorthog;
    ctrlC.progress = ctrl.progress;
    ctrlC.snapCtrl = CReflect(ctrl.snapCtrl);
    return ctrlC;
}

inline PseudospecCtrl<float> CReflect( const ElPseudospecCtrl_s& ctrlC )
{
    PseudospecCtrl<float> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}
inline PseudospecCtrl<double> CReflect( const ElPseudospecCtrl_d& ctrlC )
{
    PseudospecCtrl<double> ctrl;
    ctrl.norm = CReflect(ctrlC.norm);
    ctrl.blockWidth = ctrlC.norm;
    ctrl.schur = ctrlC.schur;
    ctrl.forceComplexSchur = ctrlC.forceComplexSchur;
    ctrl.forceComplexPs = ctrlC.forceComplexPs;
    ctrl.schurCtrl = CReflect(ctrlC.schurCtrl);
    ctrl.maxIts = ctrlC.maxIts;
    ctrl.tol = ctrlC.tol;
    ctrl.deflate = ctrlC.deflate;
    ctrl.arnoldi = ctrlC.arnoldi;
    ctrl.basisSize = ctrlC.basisSize;
    ctrl.reorthog = ctrlC.reorthog;
    ctrl.progress = ctrlC.progress;
    ctrl.snapCtrl = CReflect(ctrlC.snapCtrl);
    return ctrl;
}

inline ElSpectralBox_s CReflect( const SpectralBox<float>& box )
{
    ElSpectralBox_s boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline ElSpectralBox_d CReflect( const SpectralBox<double>& box )
{
    ElSpectralBox_d boxC;
    boxC.center = CReflect(box.center);
    boxC.realWidth = box.realWidth;
    boxC.imagWidth = box.imagWidth;
    return boxC;
}

inline SpectralBox<float> CReflect( const ElSpectralBox_s& boxC )
{
    SpectralBox<float> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

inline SpectralBox<double> CReflect( const ElSpectralBox_d& boxC )
{
    SpectralBox<double> box;
    box.center = CReflect(boxC.center);
    box.realWidth = CReflect(boxC.realWidth);
    box.imagWidth = CReflect(boxC.imagWidth);
    return box;
}

// Solvers
// ^^^^^^^
inline ElTikhonovAlg CReflect( TikhonovAlg alg )
{ return static_cast<ElTikhonovAlg>(alg); }
inline TikhonovAlg CReflect( ElTikhonovAlg alg )
{ return static_cast<TikhonovAlg>(alg); }

inline ElRidgeAlg CReflect( RidgeAlg alg )
{ return static_cast<ElRidgeAlg>(alg); }
inline RidgeAlg CReflect( ElRidgeAlg alg )
{ return static_cast<RidgeAlg>(alg); }

// Optimization
// ------------
inline ElRegularization CReflect( Regularization penalty )
{ return static_cast<ElRegularization>(penalty); }
inline Regularization CReflect( ElRegularization penalty )
{ return static_cast<Regularization>(penalty); }

inline ElKKTSystem CReflect( KKTSystem system )
{ return static_cast<ElKKTSystem>(system); }
inline KKTSystem CReflect( ElKKTSystem system )
{ return static_cast<KKTSystem>(system); }

/* Mehrotra's Predictor-Corrector IPM
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElMehrotraCtrl_s CReflect( const MehrotraCtrl<float>& ctrl )
{
    ElMehrotraCtrl_s ctrlC;
    ctrlC.primalInit    = ctrl.primalInit;
    ctrlC.dualInit      = ctrl.dualInit;
    ctrlC.minTol        = ctrl.minTol;
    ctrlC.targetTol     = ctrl.targetTol;
    ctrlC.maxIts        = ctrl.maxIts;
    ctrlC.maxStepRatio  = ctrl.maxStepRatio;
    ctrlC.system        = CReflect(ctrl.system);
    ctrlC.mehrotra      = ctrl.mehrotra;
    ctrlC.forceSameStep = ctrl.forceSameStep;
    ctrlC.solveCtrl     = CReflect(ctrl.solveCtrl);
    ctrlC.resolveReg    = ctrl.resolveReg;
    ctrlC.outerEquil    = ctrl.outerEquil;
    ctrlC.basisSize     = ctrl.basisSize;
    ctrlC.print         = ctrl.print;
    ctrlC.time          = ctrl.time;
    ctrlC.wSafeMaxNorm  = ctrl.wSafeMaxNorm;
    ctrlC.wMaxLimit     = ctrl.wMaxLimit;
    ctrlC.ruizEquilTol  = ctrl.ruizEquilTol;
    ctrlC.ruizMaxIter   = ctrl.ruizMaxIter;
    ctrlC.diagEquilTol  = ctrl.diagEquilTol;
    ctrlC.checkResiduals = ctrl.checkResiduals;
    return ctrlC;
}
inline ElMehrotraCtrl_d CReflect( const MehrotraCtrl<double>& ctrl )
{
    ElMehrotraCtrl_d ctrlC;
    ctrlC.primalInit    = ctrl.primalInit;
    ctrlC.dualInit      = ctrl.dualInit;
    ctrlC.minTol        = ctrl.minTol;
    ctrlC.targetTol     = ctrl.targetTol;
    ctrlC.maxIts        = ctrl.maxIts;
    ctrlC.maxStepRatio  = ctrl.maxStepRatio;
    ctrlC.system        = CReflect(ctrl.system);
    ctrlC.mehrotra      = ctrl.mehrotra;
    ctrlC.forceSameStep = ctrl.forceSameStep;
    ctrlC.solveCtrl     = CReflect(ctrl.solveCtrl);
    ctrlC.resolveReg    = ctrl.resolveReg;
    ctrlC.outerEquil    = ctrl.outerEquil;
    ctrlC.basisSize     = ctrl.basisSize;
    ctrlC.print         = ctrl.print;
    ctrlC.time          = ctrl.time;
    ctrlC.wSafeMaxNorm  = ctrl.wSafeMaxNorm;
    ctrlC.wMaxLimit     = ctrl.wMaxLimit;
    ctrlC.ruizEquilTol  = ctrl.ruizEquilTol;
    ctrlC.ruizMaxIter   = ctrl.ruizMaxIter;
    ctrlC.diagEquilTol  = ctrl.diagEquilTol;
    ctrlC.checkResiduals = ctrl.checkResiduals;
    return ctrlC;
}
inline MehrotraCtrl<float> CReflect( const ElMehrotraCtrl_s& ctrlC )
{
    MehrotraCtrl<float> ctrl;
    ctrl.primalInit    = ctrlC.primalInit;
    ctrl.dualInit      = ctrlC.dualInit;
    ctrl.minTol        = ctrlC.minTol;
    ctrl.targetTol     = ctrlC.targetTol;
    ctrl.maxIts        = ctrlC.maxIts;
    ctrl.maxStepRatio  = ctrlC.maxStepRatio;
    ctrl.system        = CReflect(ctrlC.system);
    ctrl.mehrotra      = ctrlC.mehrotra;
    ctrl.forceSameStep = ctrlC.forceSameStep;
    ctrl.solveCtrl     = CReflect(ctrlC.solveCtrl);
    ctrl.resolveReg    = ctrlC.resolveReg;
    ctrl.outerEquil    = ctrlC.outerEquil;
    ctrl.basisSize     = ctrlC.basisSize;
    ctrl.print         = ctrlC.print;
    ctrl.time          = ctrlC.time;
    ctrl.wSafeMaxNorm  = ctrlC.wSafeMaxNorm;
    ctrl.wMaxLimit     = ctrlC.wMaxLimit;
    ctrl.ruizEquilTol  = ctrlC.ruizEquilTol;
    ctrl.ruizMaxIter   = ctrlC.ruizMaxIter;
    ctrl.diagEquilTol  = ctrlC.diagEquilTol;
    ctrl.checkResiduals = ctrlC.checkResiduals;
    return ctrl;
}
inline MehrotraCtrl<double> CReflect( const ElMehrotraCtrl_d& ctrlC )
{
    MehrotraCtrl<double> ctrl;
    ctrl.primalInit    = ctrlC.primalInit;
    ctrl.dualInit      = ctrlC.dualInit;
    ctrl.minTol        = ctrlC.minTol;
    ctrl.targetTol     = ctrlC.targetTol;
    ctrl.maxIts        = ctrlC.maxIts;
    ctrl.maxStepRatio  = ctrlC.maxStepRatio;
    ctrl.system        = CReflect(ctrlC.system);
    ctrl.mehrotra      = ctrlC.mehrotra;
    ctrl.forceSameStep = ctrlC.forceSameStep;
    ctrl.solveCtrl     = CReflect(ctrlC.solveCtrl);
    ctrl.resolveReg    = ctrlC.resolveReg;
    ctrl.outerEquil    = ctrlC.outerEquil;
    ctrl.basisSize     = ctrlC.basisSize;
    ctrl.print         = ctrlC.print;
    ctrl.time          = ctrlC.time;
    ctrl.wSafeMaxNorm  = ctrlC.wSafeMaxNorm;
    ctrl.wMaxLimit     = ctrlC.wMaxLimit;
    ctrl.ruizEquilTol  = ctrlC.ruizEquilTol;
    ctrl.ruizMaxIter   = ctrlC.ruizMaxIter;
    ctrl.diagEquilTol  = ctrlC.diagEquilTol;
    ctrl.checkResiduals = ctrlC.checkResiduals;
    return ctrl;
}

/* Alternating Direction Method of Multipliers
   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElADMMCtrl_s CReflect( const ADMMCtrl<float>& ctrl )
{
    ElADMMCtrl_s ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline ElADMMCtrl_d CReflect( const ADMMCtrl<double>& ctrl )
{
    ElADMMCtrl_d ctrlC;
    ctrlC.rho     = ctrl.rho;
    ctrlC.alpha   = ctrl.alpha;
    ctrlC.maxIter = ctrl.maxIter;
    ctrlC.absTol  = ctrl.absTol;
    ctrlC.relTol  = ctrl.relTol;
    ctrlC.inv     = ctrl.inv;
    ctrlC.print   = ctrl.print;
    return ctrlC;
}
inline ADMMCtrl<float> CReflect( const ElADMMCtrl_s& ctrlC )
{
    ADMMCtrl<float> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}
inline ADMMCtrl<double> CReflect( const ElADMMCtrl_d& ctrlC )
{
    ADMMCtrl<double> ctrl;
    ctrl.rho     = ctrlC.rho;
    ctrl.alpha   = ctrlC.alpha;
    ctrl.maxIter = ctrlC.maxIter;
    ctrl.absTol  = ctrlC.absTol;
    ctrl.relTol  = ctrlC.relTol;
    ctrl.inv     = ctrlC.inv;
    ctrl.print   = ctrlC.print;
    return ctrl;
}

/* Linear programs
   ^^^^^^^^^^^^^^^ */
inline ElLPApproach CReflect( LPApproach approach )
{ return static_cast<ElLPApproach>(approach); }
inline LPApproach CReflect( ElLPApproach approach )
{ return static_cast<LPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElLPDirectCtrl_s CReflect( const lp::direct::Ctrl<float>& ctrl )
{
    ElLPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPDirectCtrl_d CReflect( const lp::direct::Ctrl<double>& ctrl )
{
    ElLPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.admmCtrl     = CReflect(ctrl.admmCtrl);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::direct::Ctrl<float> CReflect( const ElLPDirectCtrl_s& ctrlC )
{
    lp::direct::Ctrl<float> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::direct::Ctrl<double> CReflect( const ElLPDirectCtrl_d& ctrlC )
{
    lp::direct::Ctrl<double> ctrl(false);
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.admmCtrl     = CReflect(ctrlC.admmCtrl);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElLPAffineCtrl_s CReflect( const lp::affine::Ctrl<float>& ctrl )
{
    ElLPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElLPAffineCtrl_d CReflect( const lp::affine::Ctrl<double>& ctrl )
{
    ElLPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline lp::affine::Ctrl<float> CReflect( const ElLPAffineCtrl_s& ctrlC )
{
    lp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline lp::affine::Ctrl<double> CReflect( const ElLPAffineCtrl_d& ctrlC )
{
    lp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Quadratic programs
   ^^^^^^^^^^^^^^^^^^ */
inline ElQPApproach CReflect( QPApproach approach )
{ return static_cast<ElQPApproach>(approach); }
inline QPApproach CReflect( ElQPApproach approach )
{ return static_cast<QPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElQPDirectCtrl_s CReflect( const qp::direct::Ctrl<float>& ctrl )
{
    ElQPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPDirectCtrl_d CReflect( const qp::direct::Ctrl<double>& ctrl )
{
    ElQPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::direct::Ctrl<float> CReflect( const ElQPDirectCtrl_s& ctrlC )
{
    qp::direct::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::direct::Ctrl<double> CReflect( const ElQPDirectCtrl_d& ctrlC )
{
    qp::direct::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElQPAffineCtrl_s CReflect( const qp::affine::Ctrl<float>& ctrl )
{
    ElQPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElQPAffineCtrl_d CReflect( const qp::affine::Ctrl<double>& ctrl )
{
    ElQPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline qp::affine::Ctrl<float> CReflect( const ElQPAffineCtrl_s& ctrlC )
{
    qp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline qp::affine::Ctrl<double> CReflect( const ElQPAffineCtrl_d& ctrlC )
{
    qp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Second-order cone programs
   ^^^^^^^^^^^^^^^^^^^^^^^^^^ */
inline ElSOCPApproach CReflect( SOCPApproach approach )
{ return static_cast<ElSOCPApproach>(approach); }
inline SOCPApproach CReflect( ElSOCPApproach approach )
{ return static_cast<SOCPApproach>(approach); }

/* Direct conic form
   """"""""""""""""" */
inline ElSOCPDirectCtrl_s CReflect( const socp::direct::Ctrl<float>& ctrl )
{
    ElSOCPDirectCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElSOCPDirectCtrl_d CReflect( const socp::direct::Ctrl<double>& ctrl )
{
    ElSOCPDirectCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline socp::direct::Ctrl<float> CReflect( const ElSOCPDirectCtrl_s& ctrlC )
{
    socp::direct::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline socp::direct::Ctrl<double> CReflect( const ElSOCPDirectCtrl_d& ctrlC )
{
    socp::direct::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

/* Affine conic form
   """"""""""""""""" */
inline ElSOCPAffineCtrl_s CReflect( const socp::affine::Ctrl<float>& ctrl )
{
    ElSOCPAffineCtrl_s ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline ElSOCPAffineCtrl_d CReflect( const socp::affine::Ctrl<double>& ctrl )
{
    ElSOCPAffineCtrl_d ctrlC;
    ctrlC.approach     = CReflect(ctrl.approach);
    ctrlC.mehrotraCtrl = CReflect(ctrl.mehrotraCtrl);
    return ctrlC;
}
inline socp::affine::Ctrl<float> CReflect( const ElSOCPAffineCtrl_s& ctrlC )
{
    socp::affine::Ctrl<float> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}
inline socp::affine::Ctrl<double> CReflect( const ElSOCPAffineCtrl_d& ctrlC )
{
    socp::affine::Ctrl<double> ctrl;
    ctrl.approach     = CReflect(ctrlC.approach);
    ctrl.mehrotraCtrl = CReflect(ctrlC.mehrotraCtrl);
    return ctrl;
}

// Models
// ^^^^^^

// Basis Pursuit
// """""""""""""
inline ElBPADMMCtrl_s CReflect( const bp::ADMMCtrl<float>& ctrl )
{
    ElBPADMMCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.usePinv  = ctrl.usePinv;
    ctrlC.pinvTol  = ctrl.pinvTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElBPADMMCtrl_d CReflect( const bp::ADMMCtrl<double>& ctrl )
{
    ElBPADMMCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.usePinv  = ctrl.usePinv;
    ctrlC.pinvTol  = ctrl.pinvTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline bp::ADMMCtrl<float> CReflect( const ElBPADMMCtrl_s& ctrlC )
{
    bp::ADMMCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.usePinv  = ctrlC.usePinv;
    ctrl.pinvTol  = ctrlC.pinvTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline bp::ADMMCtrl<double> CReflect( const ElBPADMMCtrl_d& ctrlC )
{
    bp::ADMMCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.usePinv  = ctrlC.usePinv;
    ctrl.pinvTol  = ctrlC.pinvTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ElBPCtrl_s CReflect( const BPCtrl<float>& ctrl )
{
    ElBPCtrl_s ctrlC;
    ctrlC.useIPM      = ctrl.useIPM;
    ctrlC.useSOCP     = ctrl.useSOCP;
    ctrlC.admmCtrl    = CReflect(ctrl.admmCtrl);
    ctrlC.lpIPMCtrl   = CReflect(ctrl.lpIPMCtrl);
    ctrlC.socpIPMCtrl = CReflect(ctrl.socpIPMCtrl);
    return ctrlC;
}

inline ElBPCtrl_d CReflect( const BPCtrl<double>& ctrl )
{
    ElBPCtrl_d ctrlC;
    ctrlC.useIPM      = ctrl.useIPM;
    ctrlC.useSOCP     = ctrl.useSOCP;
    ctrlC.admmCtrl    = CReflect(ctrl.admmCtrl);
    ctrlC.lpIPMCtrl   = CReflect(ctrl.lpIPMCtrl);
    ctrlC.socpIPMCtrl = CReflect(ctrl.socpIPMCtrl);
    return ctrlC;
}

inline ElBPCtrl_c CReflect( const BPCtrl<Complex<float>>& ctrl )
{
    ElBPCtrl_c ctrlC;
    ctrlC.ipmCtrl = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElBPCtrl_z CReflect( const BPCtrl<Complex<double>>& ctrl )
{
    ElBPCtrl_z ctrlC;
    ctrlC.ipmCtrl = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline BPCtrl<float> CReflect( const ElBPCtrl_s& ctrlC )
{
    BPCtrl<float> ctrl(false);
    ctrl.useIPM      = ctrlC.useIPM;
    ctrl.useSOCP     = ctrlC.useSOCP;
    ctrl.admmCtrl    = CReflect(ctrlC.admmCtrl);
    ctrl.lpIPMCtrl   = CReflect(ctrlC.lpIPMCtrl);
    ctrl.socpIPMCtrl = CReflect(ctrlC.socpIPMCtrl);
    return ctrl;
}

inline BPCtrl<double> CReflect( const ElBPCtrl_d& ctrlC )
{
    BPCtrl<double> ctrl(false);
    ctrl.useIPM      = ctrlC.useIPM;
    ctrl.useSOCP     = ctrlC.useSOCP;
    ctrl.admmCtrl    = CReflect(ctrlC.admmCtrl);
    ctrl.lpIPMCtrl   = CReflect(ctrlC.lpIPMCtrl);
    ctrl.socpIPMCtrl = CReflect(ctrlC.socpIPMCtrl);
    return ctrl;
}

inline BPCtrl<Complex<float>> CReflect( const ElBPCtrl_c& ctrlC )
{
    BPCtrl<Complex<float>> ctrl;
    ctrl.ipmCtrl = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline BPCtrl<Complex<double>> CReflect( const ElBPCtrl_z& ctrlC )
{
    BPCtrl<Complex<double>> ctrl;
    ctrl.ipmCtrl = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}


// BPDN / LASSO
// """"""""""""

inline ElBPDNADMMCtrl_s CReflect( const bpdn::ADMMCtrl<float>& ctrl )
{
    ElBPDNADMMCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElBPDNADMMCtrl_d CReflect( const bpdn::ADMMCtrl<double>& ctrl )
{
    ElBPDNADMMCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline bpdn::ADMMCtrl<float> CReflect( const ElBPDNADMMCtrl_s& ctrlC )
{
    bpdn::ADMMCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline bpdn::ADMMCtrl<double> CReflect( const ElBPDNADMMCtrl_d& ctrlC )
{
    bpdn::ADMMCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ElBPDNCtrl_s CReflect( const BPDNCtrl<float>& ctrl )
{
    ElBPDNCtrl_s ctrlC;
    ctrlC.useIPM   = ctrl.useIPM;
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.ipmCtrl  = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElBPDNCtrl_d CReflect( const BPDNCtrl<double>& ctrl )
{
    ElBPDNCtrl_d ctrlC;
    ctrlC.useIPM   = ctrl.useIPM;
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.ipmCtrl  = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline BPDNCtrl<float> CReflect( const ElBPDNCtrl_s& ctrlC )
{
    BPDNCtrl<float> ctrl;
    ctrl.useIPM   = ctrlC.useIPM;
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.ipmCtrl  = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline BPDNCtrl<double> CReflect( const ElBPDNCtrl_d& ctrlC )
{
    BPDNCtrl<double> ctrl;
    ctrl.useIPM   = ctrlC.useIPM;
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.ipmCtrl  = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

// Non-negative Least Squares
// """"""""""""""""""""""""""

inline NNLSApproach CReflect( ElNNLSApproach approach ) 
{ return static_cast<NNLSApproach>(approach); }
inline ElNNLSApproach CReflect( NNLSApproach approach )
{ return static_cast<ElNNLSApproach>(approach); }

inline ElNNLSCtrl_s CReflect( const NNLSCtrl<float>& ctrl )
{
    ElNNLSCtrl_s ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.qpCtrl   = CReflect(ctrl.qpCtrl);
    ctrlC.socpCtrl = CReflect(ctrl.socpCtrl);
    return ctrlC;
}

inline ElNNLSCtrl_d CReflect( const NNLSCtrl<double>& ctrl )
{
    ElNNLSCtrl_d ctrlC;
    ctrlC.approach = CReflect(ctrl.approach);
    ctrlC.admmCtrl = CReflect(ctrl.admmCtrl);
    ctrlC.qpCtrl   = CReflect(ctrl.qpCtrl);
    ctrlC.socpCtrl = CReflect(ctrl.socpCtrl);
    return ctrlC;
}

inline NNLSCtrl<float> CReflect( const ElNNLSCtrl_s& ctrlC )
{
    NNLSCtrl<float> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.qpCtrl   = CReflect(ctrlC.qpCtrl);
    ctrl.socpCtrl = CReflect(ctrlC.socpCtrl);
    return ctrl;
}

inline NNLSCtrl<double> CReflect( const ElNNLSCtrl_d& ctrlC )
{
    NNLSCtrl<double> ctrl;
    ctrl.approach = CReflect(ctrlC.approach);
    ctrl.admmCtrl = CReflect(ctrlC.admmCtrl);
    ctrl.qpCtrl   = CReflect(ctrlC.qpCtrl);
    ctrl.socpCtrl = CReflect(ctrlC.socpCtrl);
    return ctrl;
}

// Non-negative Matrix Factorization
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline ElNMFCtrl_s CReflect( const NMFCtrl<float>& ctrl )
{
    ElNMFCtrl_s ctrlC;
    ctrlC.nnlsCtrl = CReflect(ctrl.nnlsCtrl);
    ctrlC.maxIter = ctrl.maxIter;
    return ctrlC;
}

inline ElNMFCtrl_d CReflect( const NMFCtrl<double>& ctrl )
{
    ElNMFCtrl_d ctrlC;
    ctrlC.nnlsCtrl = CReflect(ctrl.nnlsCtrl);
    ctrlC.maxIter = ctrl.maxIter;
    return ctrlC;
}

inline NMFCtrl<float> CReflect( const ElNMFCtrl_s& ctrlC )
{
    NMFCtrl<float> ctrl;
    ctrl.nnlsCtrl = CReflect(ctrlC.nnlsCtrl);
    ctrl.maxIter = ctrlC.maxIter;
    return ctrl;
}

inline NMFCtrl<double> CReflect( const ElNMFCtrl_d& ctrlC )
{
    NMFCtrl<double> ctrl;
    ctrl.nnlsCtrl = CReflect(ctrlC.nnlsCtrl);
    ctrl.maxIter = ctrlC.maxIter;
    return ctrl;
}

// Robust Principal Component Analysis
// """""""""""""""""""""""""""""""""""

inline ElRPCACtrl_s CReflect( const RPCACtrl<float>& ctrl )
{
    ElRPCACtrl_s ctrlC;
    ctrlC.useALM      = ctrl.useALM;
    ctrlC.usePivQR    = ctrl.usePivQR;
    ctrlC.progress    = ctrl.progress;
    ctrlC.numPivSteps = ctrl.numPivSteps;
    ctrlC.maxIts      = ctrl.maxIts;
    ctrlC.tau         = ctrl.tau;
    ctrlC.beta        = ctrl.beta;
    ctrlC.rho         = ctrl.rho;
    ctrlC.tol         = ctrl.tol;
    return ctrlC;
}

inline ElRPCACtrl_d CReflect( const RPCACtrl<double>& ctrl )
{
    ElRPCACtrl_d ctrlC;
    ctrlC.useALM      = ctrl.useALM;
    ctrlC.usePivQR    = ctrl.usePivQR;
    ctrlC.progress    = ctrl.progress;
    ctrlC.numPivSteps = ctrl.numPivSteps;
    ctrlC.maxIts      = ctrl.maxIts;
    ctrlC.tau         = ctrl.tau;
    ctrlC.beta        = ctrl.beta;
    ctrlC.rho         = ctrl.rho;
    ctrlC.tol         = ctrl.tol;
    return ctrlC;
}

inline RPCACtrl<float> CReflect( const ElRPCACtrl_s& ctrlC )
{
    RPCACtrl<float> ctrl;
    ctrl.useALM      = ctrlC.useALM;
    ctrl.usePivQR    = ctrlC.usePivQR;
    ctrl.progress    = ctrlC.progress;
    ctrl.numPivSteps = ctrlC.numPivSteps;
    ctrl.maxIts      = ctrlC.maxIts;
    ctrl.tau         = ctrlC.tau;
    ctrl.beta        = ctrlC.beta;
    ctrl.rho         = ctrlC.rho;
    ctrl.tol         = ctrlC.tol;
    return ctrl;
}

inline RPCACtrl<double> CReflect( const ElRPCACtrl_d& ctrlC )
{
    RPCACtrl<double> ctrl;
    ctrl.useALM      = ctrlC.useALM;
    ctrl.usePivQR    = ctrlC.usePivQR;
    ctrl.progress    = ctrlC.progress;
    ctrl.numPivSteps = ctrlC.numPivSteps;
    ctrl.maxIts      = ctrlC.maxIts;
    ctrl.tau         = ctrlC.tau;
    ctrl.beta        = ctrlC.beta;
    ctrl.rho         = ctrlC.rho;
    ctrl.tol         = ctrlC.tol;
    return ctrl;
}

// Sparse inverse covariance selection
// """""""""""""""""""""""""""""""""""

inline ElSparseInvCovCtrl_s CReflect( const SparseInvCovCtrl<float>& ctrl )
{
    ElSparseInvCovCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElSparseInvCovCtrl_d CReflect( const SparseInvCovCtrl<double>& ctrl )
{
    ElSparseInvCovCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.alpha    = ctrl.alpha;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.absTol   = ctrl.absTol;
    ctrlC.relTol   = ctrl.relTol;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline SparseInvCovCtrl<float> CReflect( const ElSparseInvCovCtrl_s& ctrlC )
{
    SparseInvCovCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline SparseInvCovCtrl<double> CReflect( const ElSparseInvCovCtrl_d& ctrlC )
{
    SparseInvCovCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.alpha    = ctrlC.alpha;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.absTol   = ctrlC.absTol;
    ctrl.relTol   = ctrlC.relTol;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* Model Fit
   """"""""" */
inline ElModelFitCtrl_s CReflect( const ModelFitCtrl<float>& ctrl )
{
    ElModelFitCtrl_s ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ElModelFitCtrl_d CReflect( const ModelFitCtrl<double>& ctrl )
{
    ElModelFitCtrl_d ctrlC;
    ctrlC.rho      = ctrl.rho;
    ctrlC.maxIter  = ctrl.maxIter;
    ctrlC.inv      = ctrl.inv;
    ctrlC.progress = ctrl.progress;
    return ctrlC;
}

inline ModelFitCtrl<float> CReflect( const ElModelFitCtrl_s& ctrlC )
{
    ModelFitCtrl<float> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

inline ModelFitCtrl<double> CReflect( const ElModelFitCtrl_d& ctrlC )
{
    ModelFitCtrl<double> ctrl;
    ctrl.rho      = ctrlC.rho;
    ctrl.maxIter  = ctrlC.maxIter;
    ctrl.inv      = ctrlC.inv;
    ctrl.progress = ctrlC.progress;
    return ctrl;
}

/* Support Vector Machine
   """""""""""""""""""""" */
inline ElSVMCtrl_s CReflect( const SVMCtrl<float>& ctrl )
{
    ElSVMCtrl_s ctrlC;
    ctrlC.useIPM       = ctrl.useIPM;
    ctrlC.modelFitCtrl = CReflect(ctrl.modelFitCtrl);
    ctrlC.ipmCtrl      = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline ElSVMCtrl_d CReflect( const SVMCtrl<double>& ctrl )
{
    ElSVMCtrl_d ctrlC;
    ctrlC.useIPM       = ctrl.useIPM;
    ctrlC.modelFitCtrl = CReflect(ctrl.modelFitCtrl);
    ctrlC.ipmCtrl      = CReflect(ctrl.ipmCtrl);
    return ctrlC;
}

inline SVMCtrl<float> CReflect( const ElSVMCtrl_s& ctrlC )
{
    SVMCtrl<float> ctrl;
    ctrl.useIPM       = ctrlC.useIPM;
    ctrl.modelFitCtrl = CReflect(ctrlC.modelFitCtrl);
    ctrl.ipmCtrl      = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

inline SVMCtrl<double> CReflect( const ElSVMCtrl_d& ctrlC )
{
    SVMCtrl<double> ctrl;
    ctrl.useIPM       = ctrlC.useIPM;
    ctrl.modelFitCtrl = CReflect(ctrlC.modelFitCtrl);
    ctrl.ipmCtrl      = CReflect(ctrlC.ipmCtrl);
    return ctrl;
}

// Lattice
// -------
inline ElLLLVariant CReflect( LLLVariant var )
{ return static_cast<ElLLLVariant>(var); }
inline LLLVariant CReflect( ElLLLVariant var )
{ return static_cast<LLLVariant>(var); }

inline LLLInfo<float> CReflect( const ElLLLInfo_s& infoC )
{
    LLLInfo<float> info;
    info.delta = infoC.delta;
    info.eta = infoC.eta;
    info.rank = infoC.rank;
    info.nullity = infoC.nullity;
    info.numSwaps = infoC.numSwaps;
    info.firstSwap = infoC.firstSwap;
    info.logVol = infoC.logVol;
    return info;
}

inline LLLInfo<double> CReflect( const ElLLLInfo_d& infoC )
{
    LLLInfo<double> info;
    info.delta = infoC.delta;
    info.eta = infoC.eta;
    info.rank = infoC.rank;
    info.nullity = infoC.nullity;
    info.numSwaps = infoC.numSwaps;
    info.firstSwap = infoC.firstSwap;
    info.logVol = infoC.logVol;
    return info;
}

inline ElLLLInfo_s CReflect( const LLLInfo<float>& info )
{
    ElLLLInfo_s infoC;
    infoC.delta = info.delta;
    infoC.eta = info.eta;
    infoC.rank = info.rank;
    infoC.nullity = info.nullity;
    infoC.numSwaps = info.numSwaps;
    infoC.firstSwap = info.firstSwap;
    infoC.logVol = info.logVol;
    return infoC;
}

inline ElLLLInfo_d CReflect( const LLLInfo<double>& info )
{
    ElLLLInfo_d infoC;
    infoC.delta = info.delta;
    infoC.eta = info.eta;
    infoC.rank = info.rank;
    infoC.nullity = info.nullity;
    infoC.numSwaps = info.numSwaps;
    infoC.firstSwap = info.firstSwap;
    infoC.logVol = info.logVol;
    return infoC;
}

inline LLLCtrl<float> CReflect( const ElLLLCtrl_s& ctrlC )
{
    LLLCtrl<float> ctrl;
    ctrl.delta = ctrlC.delta;
    ctrl.eta = ctrlC.eta;
    ctrl.variant = CReflect(ctrlC.variant);
    ctrl.recursive = ctrlC.recursive;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.presort = ctrlC.presort;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    ctrl.reorthogTol = ctrlC.reorthogTol;
    ctrl.numOrthog = ctrlC.numOrthog;
    ctrl.zeroTol = ctrlC.zeroTol;
    ctrl.blockingThresh = ctrlC.blockingThresh;
    ctrl.progress = ctrlC.progress;
    ctrl.time = ctrlC.time;
    ctrl.jumpstart = ctrlC.jumpstart;
    ctrl.startCol = ctrlC.startCol;
    return ctrl;
}

inline LLLCtrl<double> CReflect( const ElLLLCtrl_d& ctrlC )
{
    LLLCtrl<double> ctrl;
    ctrl.delta = ctrlC.delta;
    ctrl.eta = ctrlC.eta;
    ctrl.variant = CReflect(ctrlC.variant);
    ctrl.recursive = ctrlC.recursive;
    ctrl.cutoff = ctrlC.cutoff;
    ctrl.presort = ctrlC.presort;
    ctrl.smallestFirst = ctrlC.smallestFirst;
    ctrl.reorthogTol = ctrlC.reorthogTol;
    ctrl.numOrthog = ctrlC.numOrthog;
    ctrl.zeroTol = ctrlC.zeroTol;
    ctrl.blockingThresh = ctrlC.blockingThresh;
    ctrl.progress = ctrlC.progress;
    ctrl.time = ctrlC.time;
    ctrl.jumpstart = ctrlC.jumpstart;
    ctrl.startCol = ctrlC.startCol;
    return ctrl;
}

inline ElLLLCtrl_s CReflect( const LLLCtrl<float>& ctrl )
{
    ElLLLCtrl_s ctrlC;
    ctrlC.delta = ctrl.delta;
    ctrlC.eta = ctrl.eta;
    ctrlC.variant = CReflect(ctrl.variant);
    ctrlC.recursive = ctrl.recursive;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.presort = ctrl.presort;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    ctrlC.reorthogTol = ctrl.reorthogTol;
    ctrlC.numOrthog = ctrl.numOrthog;
    ctrlC.zeroTol = ctrl.zeroTol;
    ctrlC.blockingThresh = ctrl.blockingThresh;
    ctrlC.progress = ctrl.progress;
    ctrlC.time = ctrl.time;
    ctrlC.jumpstart = ctrl.jumpstart;
    ctrlC.startCol = ctrl.startCol;
    return ctrlC;
}

inline ElLLLCtrl_d CReflect( const LLLCtrl<double>& ctrl )
{
    ElLLLCtrl_d ctrlC;
    ctrlC.delta = ctrl.delta;
    ctrlC.eta = ctrl.eta;
    ctrlC.variant = CReflect(ctrl.variant);
    ctrlC.recursive = ctrl.recursive;
    ctrlC.cutoff = ctrl.cutoff;
    ctrlC.presort = ctrl.presort;
    ctrlC.smallestFirst = ctrl.smallestFirst;
    ctrlC.reorthogTol = ctrl.reorthogTol;
    ctrlC.numOrthog = ctrl.numOrthog;
    ctrlC.zeroTol = ctrl.zeroTol;
    ctrlC.blockingThresh = ctrl.blockingThresh;
    ctrlC.progress = ctrl.progress;
    ctrlC.time = ctrl.time;
    ctrlC.jumpstart = ctrl.jumpstart;
    ctrlC.startCol = ctrl.startCol;
    return ctrlC;
}

} // namespace El

#endif // ifndef EL_CREFLECT_C_HPP
