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

#define EL_RC(TYPE,INPUT) reinterpret_cast<TYPE>(INPUT)

namespace El {

template<typename T>
inline void DynamicCastCheck( T* A )
{ if( A == nullptr ) RuntimeError("Dynamic cast failed"); }

// Dist
// ----
inline Dist Reinterpret( ElDist dist )
{ return static_cast<Dist>(dist); }

inline ElDist Reinterpret( Dist dist )
{ return static_cast<ElDist>(dist); }

// Grid
// ----

inline GridOrder Reinterpret( ElGridOrderType order )
{ return static_cast<GridOrder>(order); }

inline ElGridOrderType Reinterpret( GridOrder order )
{ return static_cast<ElGridOrderType>(order); }

inline Grid* Reinterpret( ElGrid grid )
{ return EL_RC(Grid*,grid); }

inline const Grid* Reinterpret( ElConstGrid grid )
{ return EL_RC(const Grid*,grid); }

inline ElGrid Reinterpret( Grid* grid )
{ return (ElGrid)EL_RC(struct ElGrid_sDummy*,grid); }

inline ElConstGrid Reinterpret( const Grid* grid )
{ return (ElConstGrid)EL_RC(const struct ElGrid_sDummy*,grid); }

// Complex<T>
// ----------
inline Complex<float>* Reinterpret( complex_float* buffer )
{ return EL_RC(Complex<float>*,buffer); }

inline Complex<double>* Reinterpret( complex_double* buffer )
{ return EL_RC(Complex<double>*,buffer); }

inline const Complex<float>* Reinterpret( const complex_float* buffer )
{ return EL_RC(const Complex<float>*,buffer); }

inline const Complex<double>* Reinterpret( const complex_double* buffer )
{ return EL_RC(const Complex<double>*,buffer); }

// Matrix
// ------
inline Matrix<float>* Reinterpret( ElMatrix_s A )
{ return EL_RC(Matrix<float>*,A); }

inline Matrix<double>* Reinterpret( ElMatrix_d A )
{ return EL_RC(Matrix<double>*,A); }

inline Matrix<Complex<float>>* Reinterpret( ElMatrix_c A )
{ return EL_RC(Matrix<Complex<float>>*,A); }

inline Matrix<Complex<double>>* Reinterpret( ElMatrix_z A )
{ return EL_RC(Matrix<Complex<double>>*,A); }

inline const Matrix<float>* Reinterpret( ElConstMatrix_s A )
{ return EL_RC(const Matrix<float>*,A); }

inline const Matrix<double>* Reinterpret( ElConstMatrix_d A )
{ return EL_RC(const Matrix<double>*,A); }

inline const Matrix<Complex<float>>* Reinterpret( ElConstMatrix_c A )
{ return EL_RC(const Matrix<Complex<float>>*,A); }

inline const Matrix<Complex<double>>* Reinterpret( ElConstMatrix_z A )
{ return EL_RC(const Matrix<Complex<double>>*,A); }

inline ElMatrix_s Reinterpret( Matrix<float>* A )
{ return (ElMatrix_s)EL_RC(struct ElMatrix_sDummy*,A); }

inline ElMatrix_d Reinterpret( Matrix<double>* A )
{ return (ElMatrix_d)EL_RC(struct ElMatrix_dDummy*,A); }

inline ElMatrix_c Reinterpret( Matrix<Complex<float>>* A )
{ return (ElMatrix_c)EL_RC(struct ElMatrix_cDummy*,A); }

inline ElMatrix_z Reinterpret( Matrix<Complex<double>>* A )
{ return (ElMatrix_z)EL_RC(struct ElMatrix_zDummy*,A); }

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

} // namespace El

#endif // ifndef EL_REINTERPRET_C_HPP
