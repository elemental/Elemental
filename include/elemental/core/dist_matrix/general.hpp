/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DISTMATRIX_GENERAL_DECL_HPP
#define ELEM_DISTMATRIX_GENERAL_DECL_HPP

namespace elem {

template<typename T,Dist U,Dist V> 
class GeneralDistMatrix : public AbstractDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> admType;
    typedef GeneralDistMatrix<T,U,V> type;
#ifndef SWIG
    static constexpr Dist UDiag = DiagColDist<U,V>();
    static constexpr Dist VDiag = DiagRowDist<U,V>();
    static constexpr Dist UGath = GatheredDist<U>();
    static constexpr Dist VGath = GatheredDist<V>();
    static constexpr Dist UPart = PartialDist<U>();
    static constexpr Dist VPart = PartialDist<V>();
#endif

    // Constructors and destructors
    // ============================
#ifndef SWIG
    // Move constructor
    GeneralDistMatrix( type&& A );
#endif

    // Assignment and reconfiguration
    // ==============================
#ifndef SWIG
    // Move assignment
    type& operator=( type&& A );
    void RowSumScatterFrom( const DistMatrix<T,U,VGath>& A );
    void ColSumScatterFrom( const DistMatrix<T,UGath,V>& A );
    void SumScatterFrom( const DistMatrix<T,UGath,VGath>& A );
    void PartialRowSumScatterFrom( const DistMatrix<T,U,VPart>& A );
    void PartialColSumScatterFrom( const DistMatrix<T,UPart,V>& A );

    void RowSumScatterUpdate( T alpha, const DistMatrix<T,U,VGath>& A );
    void ColSumScatterUpdate( T alpha, const DistMatrix<T,UGath,V>& A );
    void SumScatterUpdate( T alpha, const DistMatrix<T,UGath,VGath>& A );
    void PartialRowSumScatterUpdate( T alpha, const DistMatrix<T,U,VPart>& A );
    void PartialColSumScatterUpdate( T alpha, const DistMatrix<T,UPart,V>& A );
#endif

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith( const elem::DistData& d, Int offset=0 ) const;
    Int DiagonalRoot( Int offset=0 ) const;
    Int DiagonalAlign( Int offset=0 ) const;
#ifndef SWIG
    void GetDiagonal( DistMatrix<T,UDiag,VDiag>& d, Int offset=0 ) const;
    void GetRealPartOfDiagonal
    ( DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 ) const;
    void GetImagPartOfDiagonal
    ( DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 ) const;

    DistMatrix<T,UDiag,VDiag> GetDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),UDiag,VDiag> GetRealPartOfDiagonal( Int offset=0 ) const;
    DistMatrix<BASE(T),UDiag,VDiag> GetImagPartOfDiagonal( Int offset=0 ) const;

    void SetDiagonal( const DistMatrix<T,UDiag,VDiag>& d, Int offset=0 );
    void SetRealPartOfDiagonal
    ( const DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 );
    void SetImagPartOfDiagonal
    ( const DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 );

    void UpdateDiagonal
    ( T alpha, const DistMatrix<T,UDiag,VDiag>& d, Int offset=0 );
    void UpdateRealPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 );
    void UpdateImagPartOfDiagonal
    ( BASE(T) alpha, const DistMatrix<BASE(T),UDiag,VDiag>& d, Int offset=0 );
#endif // ifndef SWIG

protected:
    // Construct using a particular process grid
    // =========================================
    GeneralDistMatrix( const elem::Grid& g );

    // Redistribution helper routines
    // ==============================
#ifndef SWIG
    void AllGather( DistMatrix<T,UGath,VGath>& A ) const;
    void ColAllGather( DistMatrix<T,UGath,V>& A ) const;
    void RowAllGather( DistMatrix<T,U,VGath>& A ) const;
    void PartialColAllGather( DistMatrix<T,UPart,V>& A ) const;
    void PartialRowAllGather( DistMatrix<T,U,VPart>& A ) const;

    void FilterFrom( const DistMatrix<T,UGath,VGath>& A );
    void ColFilterFrom( const DistMatrix<T,UGath,V>& A );
    void RowFilterFrom( const DistMatrix<T,U,VGath>& A );
    void PartialColFilterFrom( const DistMatrix<T,UPart,V>& A );
    void PartialRowFilterFrom( const DistMatrix<T,U,VPart>& A );
#endif

    // Diagonal helper routines
    // ========================
#ifndef SWIG
    template<typename S,class Function>
    void GetDiagonalHelper
    ( DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const DistMatrix<S,UDiag,VDiag>& d, Int offset, Function func );
#endif // ifndef SWIG

    // Friend declarations
    // ===================
#ifndef SWIG
    template<typename S,Dist J,Dist K> friend class DistMatrix;
#endif 
};

} // namespace elem

#endif // ifndef ELEM_DISTMATRIX_GENERAL_DECL_HPP
