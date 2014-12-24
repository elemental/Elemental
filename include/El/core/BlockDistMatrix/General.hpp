/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP
#define EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP

namespace El {

template<typename T,Dist U,Dist V> 
class GeneralBlockDistMatrix : public AbstractBlockDistMatrix<T>
{
public:
    // Typedefs
    // ========
    typedef AbstractBlockDistMatrix<T> absType;
    typedef GeneralBlockDistMatrix<T,U,V> type;
    static constexpr Dist UGath = Collect<U>();
    static constexpr Dist VGath = Collect<V>();
    static constexpr Dist UPart = Partial<U>();
    static constexpr Dist VPart = Partial<V>();

    static constexpr Dist UDiag = DiagCol<U,V>();
    static constexpr Dist VDiag = DiagRow<U,V>();
    static constexpr Dist UScat = PartialUnionCol<U,V>();
    static constexpr Dist VScat = PartialUnionRow<U,V>();

    // Constructors and destructors
    // ============================
    // Move constructor
    GeneralBlockDistMatrix( type&& A ) EL_NOEXCEPT;

    // Assignment and reconfiguration
    // ==============================
    // Move assignment
    type& operator=( type&& A );

    // Diagonal manipulation
    // =====================
    void GetDiagonal
    ( AbstractBlockDistMatrix<T>& d, Int offset=0 ) const override;
    void GetRealPartOfDiagonal
    ( AbstractBlockDistMatrix<Base<T>>& d, Int offset=0 ) const override;
    void GetImagPartOfDiagonal
    ( AbstractBlockDistMatrix<Base<T>>& d, Int offset=0 ) const override;

    BlockDistMatrix<T,UDiag,VDiag> GetDiagonal( Int offset=0 ) const;
    BlockDistMatrix<Base<T>,UDiag,VDiag> 
    GetRealPartOfDiagonal( Int offset=0 ) const;
    BlockDistMatrix<Base<T>,UDiag,VDiag> 
    GetImagPartOfDiagonal( Int offset=0 ) const;

    void SetDiagonal
    ( const AbstractBlockDistMatrix<T>& d, Int offset=0 ) override;
    void SetRealPartOfDiagonal
    ( const AbstractBlockDistMatrix<Base<T>>& d, Int offset=0 ) override;
    void SetImagPartOfDiagonal
    ( const AbstractBlockDistMatrix<Base<T>>& d, Int offset=0 ) override;

    void UpdateDiagonal
    ( T alpha, const AbstractBlockDistMatrix<T>& d, 
      Int offset=0 ) override;
    void UpdateRealPartOfDiagonal
    ( Base<T> alpha, const AbstractBlockDistMatrix<Base<T>>& d, 
      Int offset=0 ) override;
    void UpdateImagPartOfDiagonal
    ( Base<T> alpha, const AbstractBlockDistMatrix<Base<T>>& d, 
      Int offset=0 ) override;

protected:

    // Private constructors
    // ====================

    // Inherited constructors are part of C++11 but not yet widely supported.
    //using AbstractDistMatrix<T>::AbstractDistMatrix;

    GeneralBlockDistMatrix( const El::Grid& g=DefaultGrid(), Int root=0 );
    GeneralBlockDistMatrix
    ( const El::Grid& g, Int blockHeight, Int blockWidth, Int root=0 );

    // Diagonal helper routines
    // ========================
    template<typename S,class Function>
    void GetDiagonalHelper
    ( AbstractBlockDistMatrix<S>& d, Int offset, Function func ) const;
    template<typename S,class Function>
    void SetDiagonalHelper
    ( const AbstractBlockDistMatrix<S>& d, Int offset, Function func );

    // Friend declarations
    // ===================
    template<typename S,Dist J,Dist K> friend class DistMatrix;
    template<typename S,Dist J,Dist K> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_BLOCKDISTMATRIX_GENERAL_DECL_HPP
