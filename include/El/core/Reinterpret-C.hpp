/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_REINTERPRET_C_HPP
#define EL_REINTERPRET_C_HPP

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
struct CReflect { typedef T type; };

// ElInt and Int are typedef's
/*
template<> struct CReflect<ElInt> { typedef Int type; };
template<> struct CReflect<Int> { typedef ElInt type; };
*/

template<> struct CReflect<complex_float> { typedef Complex<float> type; };
template<> struct CReflect<complex_double> { typedef Complex<double> type; };

template<> struct CReflect<Complex<float>> { typedef complex_float type; };
template<> struct CReflect<Complex<double>> { typedef complex_double type; };

#define CREFLECT(T) typename CReflect<T>::type

template<typename T>
inline void DynamicCastCheck( T* A )
{ if( A == nullptr ) RuntimeError("Dynamic cast failed"); }

inline Orientation Reinterpret( ElOrientation orient ) 
{ return static_cast<Orientation>(orient); }
inline ElOrientation Reinterpret( Orientation orient )
{ return static_cast<ElOrientation>(orient); }

inline LeftOrRight Reinterpret( ElLeftOrRight side )
{ return static_cast<LeftOrRight>(side); }
inline ElLeftOrRight Reinterpret( LeftOrRight side )
{ return static_cast<ElLeftOrRight>(side); }

inline UpperOrLower Reinterpret( ElUpperOrLower uplo )
{ return static_cast<UpperOrLower>(uplo); }
inline ElUpperOrLower Reinterpret( UpperOrLower uplo )
{ return static_cast<ElUpperOrLower>(uplo); }

inline UnitOrNonUnit Reinterpret( ElUnitOrNonUnit diag )
{ return static_cast<UnitOrNonUnit>(diag); }
inline ElUnitOrNonUnit Reinterpret( UnitOrNonUnit diag )
{ return static_cast<ElUnitOrNonUnit>(diag); }

// Dist
// ----
inline Dist   Reinterpret( ElDist dist ) { return static_cast<  Dist>(dist); }
inline ElDist Reinterpret(   Dist dist ) { return static_cast<ElDist>(dist); }

// Grid
// ----
inline   GridOrder     Reinterpret( ElGridOrderType order )
{ return static_cast<  GridOrder    >(order); }
inline ElGridOrderType Reinterpret(   GridOrder     order )
{ return static_cast<ElGridOrderType>(order); }

inline Grid* Reinterpret( ElGrid grid )
{ return EL_RC(Grid*,grid); }
inline ElGrid Reinterpret( Grid* grid )
{ return (ElGrid)EL_RC(struct ElGrid_sDummy*,grid); }

inline const Grid* Reinterpret( ElConstGrid grid )
{ return EL_RC(const Grid*,grid); }
inline ElConstGrid Reinterpret( const Grid* grid )
{ return (ElConstGrid)EL_RC(const struct ElGrid_sDummy*,grid); }

// Complex<T>
// ----------
inline complex_float* Reinterpret( Complex<float>* buffer )
{ return EL_RC(complex_float*,buffer); }

inline complex_double* Reinterpret( Complex<double>* buffer )
{ return EL_RC(complex_double*,buffer); }

inline const complex_float* Reinterpret( const Complex<float>* buffer )
{ return EL_RC(const complex_float*,buffer); }

inline const complex_double* Reinterpret( const Complex<double>* buffer )
{ return EL_RC(const complex_double*,buffer); }

inline Complex<float>* Reinterpret( complex_float* buffer )
{ return EL_RC(Complex<float>*,buffer); }

inline Complex<double>* Reinterpret( complex_double* buffer )
{ return EL_RC(Complex<double>*,buffer); }

inline const Complex<float>* Reinterpret( const complex_float* buffer )
{ return EL_RC(const Complex<float>*,buffer); }

inline const Complex<double>* Reinterpret( const complex_double* buffer )
{ return EL_RC(const Complex<double>*,buffer); }

inline Complex<float> Reinterpret( complex_float alpha )
{ return Complex<float>(alpha.real,alpha.imag); }

inline Complex<double> Reinterpret( complex_double alpha )
{ return Complex<double>(alpha.real,alpha.imag); }

inline complex_float Reinterpret( Complex<float> alpha )
{ complex_float beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

inline complex_double Reinterpret( Complex<double> alpha )
{ complex_double beta; beta.real = alpha.real(); beta.imag = alpha.imag();
  return beta; }

// Analogues for real variables and integers
// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
inline Int Reinterpret( Int alpha ) { return alpha; }
inline Int* Reinterpret( Int* buffer ) { return buffer; }
inline const Int* Reinterpret( const Int* buffer ) { return buffer; }
/*
inline ElInt Reinterpret( Int alpha ) { return alpha; }

inline ElInt* Reinterpret( Int*   buffer ) { return buffer; }
inline Int*   Reinterpret( ElInt* buffer ) { return buffer; }

inline const ElInt* Reinterpret( const Int*   buffer ) { return buffer; }
inline const Int*   Reinterpret( const ElInt* buffer ) { return buffer; }
*/

inline float Reinterpret( float alpha) { return alpha; }
inline double Reinterpret( double alpha ) { return alpha; }

inline float* Reinterpret( float* buffer ) { return buffer; }
inline double* Reinterpret( double* buffer ) { return buffer; }

inline const float* Reinterpret( const float* buffer ) { return buffer; }
inline const double* Reinterpret( const double* buffer ) { return buffer; }

inline ValueInt<Int> Reinterpret( ElValueInt_i entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_i Reinterpret( ValueInt<Int> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<float> Reinterpret( ElValueInt_s entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_s Reinterpret( ValueInt<float> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<double> Reinterpret( ElValueInt_d entryC )
{ return {entryC.value,entryC.index}; }
inline ElValueInt_d Reinterpret( ValueInt<double> entry )
{ return {entry.value,entry.index}; }

inline ValueInt<Complex<float>> Reinterpret( ElValueInt_c entryC )
{ return {Reinterpret(entryC.value),entryC.index}; }
inline ElValueInt_c Reinterpret( ValueInt<Complex<float>> entry )
{ return {Reinterpret(entry.value),entry.index}; }

inline ValueInt<Complex<double>> Reinterpret( ElValueInt_z entryC )
{ return {Reinterpret(entryC.value),entryC.index}; }
inline ElValueInt_z Reinterpret( ValueInt<Complex<double>> entry )
{ return {Reinterpret(entry.value),entry.index}; }

inline ValueIntPair<Int> Reinterpret( ElValueIntPair_i entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_i Reinterpret( ValueIntPair<Int> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<float> Reinterpret( ElValueIntPair_s entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_s Reinterpret( ValueIntPair<float> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<double> Reinterpret( ElValueIntPair_d entryC )
{ return {entryC.value,{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_d Reinterpret( ValueIntPair<double> entry )
{ return {entry.value,{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<Complex<float>> Reinterpret( ElValueIntPair_c entryC )
{ return {Reinterpret(entryC.value),{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_c Reinterpret( ValueIntPair<Complex<float>> entry )
{ return {Reinterpret(entry.value),{entry.indices[0],entry.indices[1]}}; }

inline ValueIntPair<Complex<double>> Reinterpret( ElValueIntPair_z entryC )
{ return {Reinterpret(entryC.value),{entryC.indices[0],entryC.indices[1]}}; }
inline ElValueIntPair_z Reinterpret( ValueIntPair<Complex<double>> entry )
{ return {Reinterpret(entry.value),{entry.indices[0],entry.indices[1]}}; }

// Matrix
// ------
inline Matrix<Int>* Reinterpret( ElMatrix_i A )
{ return EL_RC(Matrix<Int>*,A); }

inline Matrix<float>* Reinterpret( ElMatrix_s A )
{ return EL_RC(Matrix<float>*,A); }

inline Matrix<double>* Reinterpret( ElMatrix_d A )
{ return EL_RC(Matrix<double>*,A); }

inline Matrix<Complex<float>>* Reinterpret( ElMatrix_c A )
{ return EL_RC(Matrix<Complex<float>>*,A); }

inline Matrix<Complex<double>>* Reinterpret( ElMatrix_z A )
{ return EL_RC(Matrix<Complex<double>>*,A); }

inline const Matrix<Int>* Reinterpret( ElConstMatrix_i A )
{ return EL_RC(const Matrix<Int>*,A); }

inline const Matrix<float>* Reinterpret( ElConstMatrix_s A )
{ return EL_RC(const Matrix<float>*,A); }

inline const Matrix<double>* Reinterpret( ElConstMatrix_d A )
{ return EL_RC(const Matrix<double>*,A); }

inline const Matrix<Complex<float>>* Reinterpret( ElConstMatrix_c A )
{ return EL_RC(const Matrix<Complex<float>>*,A); }

inline const Matrix<Complex<double>>* Reinterpret( ElConstMatrix_z A )
{ return EL_RC(const Matrix<Complex<double>>*,A); }

inline ElMatrix_i Reinterpret( Matrix<Int>* A )
{ return (ElMatrix_i)EL_RC(struct ElMatrix_iDummy*,A); }

inline ElMatrix_s Reinterpret( Matrix<float>* A )
{ return (ElMatrix_s)EL_RC(struct ElMatrix_sDummy*,A); }

inline ElMatrix_d Reinterpret( Matrix<double>* A )
{ return (ElMatrix_d)EL_RC(struct ElMatrix_dDummy*,A); }

inline ElMatrix_c Reinterpret( Matrix<Complex<float>>* A )
{ return (ElMatrix_c)EL_RC(struct ElMatrix_cDummy*,A); }

inline ElMatrix_z Reinterpret( Matrix<Complex<double>>* A )
{ return (ElMatrix_z)EL_RC(struct ElMatrix_zDummy*,A); }

inline ElConstMatrix_i Reinterpret( const Matrix<Int>* A )
{ return (ElConstMatrix_i)EL_RC(const struct ElMatrix_iDummy*,A); }

inline ElConstMatrix_s Reinterpret( const Matrix<float>* A )
{ return (ElConstMatrix_s)EL_RC(const struct ElMatrix_sDummy*,A); }

inline ElConstMatrix_d Reinterpret( const Matrix<double>* A )
{ return (ElConstMatrix_d)EL_RC(const struct ElMatrix_dDummy*,A); }

inline ElConstMatrix_c Reinterpret( const Matrix<Complex<float>>* A )
{ return (ElConstMatrix_c)EL_RC(const struct ElMatrix_cDummy*,A); }

inline ElConstMatrix_z Reinterpret( const Matrix<Complex<double>>* A )
{ return (ElConstMatrix_z)EL_RC(const struct ElMatrix_zDummy*,A); }

// AbstractDistMatrix
// ------------------
inline AbstractDistMatrix<Int>* 
Reinterpret( ElDistMatrix_i A )
{ return EL_RC(AbstractDistMatrix<Int>*,A); }

inline AbstractDistMatrix<float>* 
Reinterpret( ElDistMatrix_s A )
{ return EL_RC(AbstractDistMatrix<float>*,A); }

inline AbstractDistMatrix<double>* 
Reinterpret( ElDistMatrix_d A )
{ return EL_RC(AbstractDistMatrix<double>*,A); }

inline AbstractDistMatrix<Complex<float>>* 
Reinterpret( ElDistMatrix_c A )
{ return EL_RC(AbstractDistMatrix<Complex<float>>*,A); }

inline AbstractDistMatrix<Complex<double>>* 
Reinterpret( ElDistMatrix_z A )
{ return EL_RC(AbstractDistMatrix<Complex<double>>*,A); }

inline const AbstractDistMatrix<Int>* 
Reinterpret( ElConstDistMatrix_i A )
{ return EL_RC(const AbstractDistMatrix<Int>*,A); }

inline const AbstractDistMatrix<float>* 
Reinterpret( ElConstDistMatrix_s A )
{ return EL_RC(const AbstractDistMatrix<float>*,A); }

inline const AbstractDistMatrix<double>* 
Reinterpret( ElConstDistMatrix_d A )
{ return EL_RC(const AbstractDistMatrix<double>*,A); }

inline const AbstractDistMatrix<Complex<float>>* 
Reinterpret( ElConstDistMatrix_c A )
{ return EL_RC(const AbstractDistMatrix<Complex<float>>*,A); }

inline const AbstractDistMatrix<Complex<double>>* 
Reinterpret( ElConstDistMatrix_z A )
{ return EL_RC(const AbstractDistMatrix<Complex<double>>*,A); }

inline ElDistMatrix_i
Reinterpret( AbstractDistMatrix<Int>* A )
{ return (ElDistMatrix_i)EL_RC(struct ElDistMatrix_iDummy*,A); }

inline ElDistMatrix_s 
Reinterpret( AbstractDistMatrix<float>* A )
{ return (ElDistMatrix_s)EL_RC(struct ElDistMatrix_sDummy*,A); }

inline ElDistMatrix_d 
Reinterpret( AbstractDistMatrix<double>* A )
{ return (ElDistMatrix_d)EL_RC(struct ElDistMatrix_dDummy*,A); }

inline ElDistMatrix_c 
Reinterpret( AbstractDistMatrix<Complex<float>>* A )
{ return (ElDistMatrix_c)EL_RC(struct ElDistMatrix_cDummy*,A); }

inline ElDistMatrix_z 
Reinterpret( AbstractDistMatrix<Complex<double>>* A )
{ return (ElDistMatrix_z)EL_RC(struct ElDistMatrix_zDummy*,A); }

inline ElConstDistMatrix_i
Reinterpret( const AbstractDistMatrix<Int>* A )
{ return (ElConstDistMatrix_i)EL_RC(const struct ElDistMatrix_iDummy*,A); }

inline ElConstDistMatrix_s 
Reinterpret( const AbstractDistMatrix<float>* A )
{ return (ElConstDistMatrix_s)EL_RC(const struct ElDistMatrix_sDummy*,A); }

inline ElConstDistMatrix_d 
Reinterpret( const AbstractDistMatrix<double>* A )
{ return (ElConstDistMatrix_d)EL_RC(const struct ElDistMatrix_dDummy*,A); }

inline ElConstDistMatrix_c 
Reinterpret( const AbstractDistMatrix<Complex<float>>* A )
{ return (ElConstDistMatrix_c)EL_RC(const struct ElDistMatrix_cDummy*,A); }

inline ElConstDistMatrix_z 
Reinterpret( const AbstractDistMatrix<Complex<double>>* A )
{ return (ElConstDistMatrix_z)EL_RC(const struct ElDistMatrix_zDummy*,A); }

inline ElDistData Reinterpret( const DistData& data )
{
    ElDistData distData;
    distData.colDist = Reinterpret(data.colDist);
    distData.rowDist = Reinterpret(data.rowDist);
    distData.colAlign = data.colAlign;
    distData.rowAlign = data.rowAlign;
    distData.root = data.root;
    distData.grid = Reinterpret(data.grid);
    return distData;
}

inline DistData Reinterpret( const ElDistData distData )
{
    DistData data;
    data.colDist = Reinterpret(distData.colDist);
    data.rowDist = Reinterpret(distData.rowDist);
    data.colAlign = distData.colAlign;
    data.rowAlign = distData.rowAlign;
    data.root = distData.root;
    data.grid = Reinterpret(distData.grid);
    return data;
}

// BLAS-like
// ---------
inline ElGemmAlgorithm Reinterpret( GemmAlgorithm alg )
{ return static_cast<ElGemmAlgorithm>(alg); }
inline GemmAlgorithm Reinterpret( ElGemmAlgorithm alg )
{ return static_cast<GemmAlgorithm>(alg); }

// LAPACK-like
// -----------

// Condensed form
// ^^^^^^^^^^^^^^
inline ElHermitianTridiagApproach 
Reinterpret( HermitianTridiagApproach approach )
{ return static_cast<ElHermitianTridiagApproach>( approach ); }

inline HermitianTridiagApproach 
Reinterpret( ElHermitianTridiagApproach approach )
{ return static_cast<HermitianTridiagApproach>( approach ); }

inline ElHermitianTridiagCtrl
Reinterpret( HermitianTridiagCtrl ctrl )
{ 
    ElHermitianTridiagCtrl ctrlC;
    ctrlC.approach = Reinterpret(ctrl.approach);
    ctrlC.order = Reinterpret(ctrl.order);
    return ctrlC;
}

inline HermitianTridiagCtrl
Reinterpret( ElHermitianTridiagCtrl ctrlC )
{ 
    HermitianTridiagCtrl ctrl;
    ctrl.approach = Reinterpret(ctrlC.approach);
    ctrl.order = Reinterpret(ctrlC.order);
    return ctrl;
}

// Factorizations
// ^^^^^^^^^^^^^^
inline ElLDLPivotType Reinterpret( LDLPivotType pivotType )
{ return static_cast<ElLDLPivotType>( pivotType ); }

inline LDLPivotType Reinterpret( ElLDLPivotType pivotType )
{ return static_cast<LDLPivotType>( pivotType ); }

inline ElLDLPivot Reinterpret( LDLPivot pivot )
{
    ElLDLPivot pivotC;
    pivotC.nb = pivot.nb;
    pivotC.from[0] = pivot.from[0];
    pivotC.from[1] = pivot.from[1];
    return pivotC;
}

inline LDLPivot Reinterpret( ElLDLPivot pivotC )
{
    LDLPivot pivot;
    pivot.nb = pivotC.nb;
    pivot.from[0] = pivotC.from[0];
    pivot.from[1] = pivotC.from[1];
    return pivot;
}

inline ElInertiaType Reinterpret( InertiaType inertia )
{ 
    ElInertiaType inertiaC;
    inertiaC.numPositive = inertia.numPositive;
    inertiaC.numNegative = inertia.numNegative;
    inertiaC.numZero = inertia.numZero;
    return inertiaC;
}

inline InertiaType Reinterpret( ElInertiaType inertiaC )
{ 
    InertiaType inertia;
    inertia.numPositive = inertiaC.numPositive;
    inertia.numNegative = inertiaC.numNegative;
    inertia.numZero = inertiaC.numZero;
    return inertia;
}

inline ElQRCtrl_s Reinterpret( QRCtrl<float> ctrl )
{ 
    ElQRCtrl_s ctrlC;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    return ctrlC;
}
inline ElQRCtrl_d Reinterpret( QRCtrl<double> ctrl )
{ 
    ElQRCtrl_d ctrlC;
    ctrlC.boundRank = ctrl.boundRank;
    ctrlC.maxRank = ctrl.maxRank;
    ctrlC.adaptive = ctrl.adaptive;
    ctrlC.tol = ctrl.tol;
    ctrlC.alwaysRecomputeNorms = ctrl.alwaysRecomputeNorms;
    return ctrlC;
}

inline QRCtrl<float> Reinterpret( ElQRCtrl_s ctrlC )
{ 
    QRCtrl<float> ctrl;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    return ctrl;
}
inline QRCtrl<double> Reinterpret( ElQRCtrl_d ctrlC )
{ 
    QRCtrl<double> ctrl;
    ctrl.boundRank = ctrlC.boundRank;
    ctrl.maxRank = ctrlC.maxRank;
    ctrl.adaptive = ctrlC.adaptive;
    ctrl.tol = ctrlC.tol;
    ctrl.alwaysRecomputeNorms = ctrlC.alwaysRecomputeNorms;
    return ctrl;
}

} // namespace El

#endif // ifndef EL_REINTERPRET_C_HPP
