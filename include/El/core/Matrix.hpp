/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRIX_HPP
#define EL_MATRIX_HPP

namespace El {

template<typename Ring>
struct OrientedMatrix
{
    Orientation orient;
    const Matrix<Ring>& matrix;

    OrientedMatrix( const Matrix<Ring>& A, Orientation orientation=NORMAL )
    : orient(orientation), matrix(A)
    { }
};

template<typename Ring>
class Matrix
{
public:    
    // Constructors and destructors
    // ============================
    // Create a 0x0 matrix
    Matrix( bool fixed=false );
    // Create a matrix with the specified dimensions
    Matrix( Int height, Int width, bool fixed=false );
    // Create a matrix with the specified dimensions and leading dimension
    Matrix( Int height, Int width, Int ldim, bool fixed=false );
    // Construct a matrix around an existing (possibly immutable) buffer
    Matrix
    ( Int height, Int width, const Ring* buffer, Int ldim, bool fixed=false );
    Matrix( Int height, Int width, Ring* buffer, Int ldim, bool fixed=false );
    // Create a copy of a matrix
    Matrix( const Matrix<Ring>& A );
    // Move the metadata from a given matrix
    Matrix( Matrix<Ring>&& A ) EL_NOEXCEPT;
    // Destructor
    ~Matrix();

    // Assignment and reconfiguration
    // ==============================

    void Empty();
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    // Reconfigure around the given buffer, but do not assume ownership
    void Attach( Int height, Int width, Ring* buffer, Int ldim );
    void LockedAttach( Int height, Int width, const Ring* buffer, Int ldim );
    // Reconfigure around the given buffer and assume ownership
    void Control( Int height, Int width, Ring* buffer, Int ldim );

    // Operator overloading
    // ====================

    // Return a view
    // -------------
          Matrix<Ring> operator()( Range<Int> I, Range<Int> J );
    const Matrix<Ring> operator()( Range<Int> I, Range<Int> J ) const;

    // Make a copy
    // -----------
    const Matrix<Ring>& operator=( const Matrix<Ring>& A );

    // Move assignment
    // ---------------
    Matrix<Ring>& operator=( Matrix<Ring>&& A );

    // Rescaling
    // ---------
    const Matrix<Ring>& operator*=( Ring alpha );

    // Addition/substraction
    // ---------------------
    const Matrix<Ring>& operator+=( const Matrix<Ring>& A );
    const Matrix<Ring>& operator-=( const Matrix<Ring>& A );

    // Basic queries
    // =============
    Int Height() const;
    Int Width() const;
    Int LDim() const;
    Int MemorySize() const;
    Int DiagonalLength( Int offset=0 ) const;
    Ring* Buffer();
    const Ring* LockedBuffer() const;
    Ring* Buffer( Int i, Int j );
    const Ring* LockedBuffer( Int i, Int j ) const;
    bool Viewing()   const;
    bool FixedSize() const;
    bool Locked()    const;

    // Single-entry manipulation
    // =========================
    Ring Get( Int i, Int j ) const;
    Base<Ring> GetRealPart( Int i, Int j ) const;
    Base<Ring> GetImagPart( Int i, Int j ) const;
    void Set( Int i, Int j, Ring alpha );
    void Set( const Entry<Ring>& entry );
    void SetRealPart( Int i, Int j, Base<Ring> alpha );
    void SetImagPart( Int i, Int j, Base<Ring> alpha );
    void SetRealPart( const Entry<Base<Ring>>& entry );
    void SetImagPart( const Entry<Base<Ring>>& entry );
    void Update( Int i, Int j, Ring alpha );
    void Update( const Entry<Ring>& entry );
    void UpdateRealPart( Int i, Int j, Base<Ring> alpha );
    void UpdateImagPart( Int i, Int j, Base<Ring> alpha );
    void UpdateRealPart( const Entry<Base<Ring>>& entry );
    void UpdateImagPart( const Entry<Base<Ring>>& entry );
    void MakeReal( Int i, Int j );
    void Conjugate( Int i, Int j );

    // Orientation
    // ===========
    OrientedMatrix<Ring> N() const
    { return OrientedMatrix<Ring>(*this,NORMAL); }

    OrientedMatrix<Ring> T() const
    { return OrientedMatrix<Ring>(*this,TRANSPOSE); }

    OrientedMatrix<Ring> H() const
    { return OrientedMatrix<Ring>(*this,ADJOINT); }

    OrientedMatrix<Ring> Orient( Orientation orient ) const
    { return OrientedMatrix<Ring>(*this,orient); }

private:
    // Member variables
    // ================
    ViewType viewType_;
    Int height_, width_, ldim_;
    const Ring* data_;
    Memory<Ring> memory_;

    // Exchange metadata with another matrix
    // =====================================
    void ShallowSwap( Matrix<Ring>& A );

    // Reconfigure without error-checking
    // ==================================
    void Empty_();
    void Resize_( Int height, Int width );
    void Resize_( Int height, Int width, Int ldim );
    void Control_( Int height, Int width, Ring* buffer, Int ldim );
    void Attach_( Int height, Int width, Ring* buffer, Int ldim );
    void LockedAttach_( Int height, Int width, const Ring* buffer, Int ldim );

    // Return a reference to a single entry without error-checking
    // ===========================================================
    const Ring& Get_( Int i, Int j ) const;
    Ring& Set_( Int i, Int j );

    // Assertions
    // ==========
    void ComplainIfReal() const;
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions( Int height, Int width, Int ldim ) const;
    void AssertValidEntry( Int i, Int j ) const;
   
    // Friend declarations
    // ===================
    template <typename Ring2> friend class Matrix;
    template <typename Ring2> friend class AbstractDistMatrix;
    template <typename Ring2> friend class AbstractBlockDistMatrix;
    template <typename Ring2,Dist U,Dist V> friend class GeneralDistMatrix;
    template <typename Ring2,Dist U,Dist V> friend class GeneralBlockDistMatrix;
    template <typename Ring2,Dist U,Dist V> friend class DistMatrix;
    template <typename Ring2,Dist U,Dist V> friend class BlockDistMatrix;
};

} // namespace El

#endif // ifndef EL_MATRIX_HPP
