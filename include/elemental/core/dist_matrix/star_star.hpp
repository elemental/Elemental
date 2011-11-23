/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP
#define ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP 1

namespace elemental {

// Partial specialization to A[* ,* ].
//
// The entire matrix is replicated across all processes.
template<typename T>
class DistMatrix<T,STAR,STAR> : public AbstractDistMatrix<T>
{
public:
    // Create a 0 x 0 distributed matrix
    DistMatrix( const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix
    DistMatrix( int height, int width, const elemental::Grid& g=DefaultGrid() );

    // Create a height x width distributed matrix with specified alignments
    // and leading dimension
    DistMatrix( int height, int width, int ldim, const elemental::Grid& g );

    // View a constant distributed matrix's buffer
    DistMatrix
    ( int height, int width, const T* buffer, int ldim, 
      const elemental::Grid& g );

    // View a mutable distributed matrix's buffer
    DistMatrix
    ( int height, int width, T* buffer, int ldim, const elemental::Grid& g );

    // Create a copy of distributed matrix A
    template<Distribution U,Distribution V>
    DistMatrix( const DistMatrix<T,U,V>& A );

    ~DistMatrix();

    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MC,MR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MC,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MD,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MD>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MR,MC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,MR,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,MC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,VC,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,VC>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,VR,STAR>& A );
    const DistMatrix<T,STAR,STAR>& operator=( const DistMatrix<T,STAR,VR>& A );
    const DistMatrix<T,STAR,STAR>& 
    operator=( const DistMatrix<T,STAR,STAR>& A );

    //------------------------------------------------------------------------//
    // Fulfillments of abstract virtual func's from AbstractDistMatrix        //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //

    virtual void SetGrid( const elemental::Grid& grid );

    virtual T Get( int i, int j ) const;
    virtual void Set( int i, int j, T alpha );
    virtual void Update( int i, int j, T alpha );

    virtual void MakeTrapezoidal
    ( Side side, Shape shape, int offset = 0 );

    virtual void ScaleTrapezoidal
    ( T alpha, Side side, Shape shape, int offset = 0 );

    virtual void ResizeTo( int height, int width );
    virtual void SetToIdentity();
    virtual void SetToRandom();
    virtual void SetToRandomHermitian();
    virtual void SetToRandomHPD();

    //
    // Routines that are only valid for complex datatypes
    //

    virtual typename RealBase<T>::type GetReal( int i, int j ) const;
    virtual typename RealBase<T>::type GetImag( int i, int j ) const;
    virtual void SetReal( int i, int j, typename RealBase<T>::type u );
    virtual void SetImag( int i, int j, typename RealBase<T>::type u );
    virtual void UpdateReal( int i, int j, typename RealBase<T>::type u );
    virtual void UpdateImag( int i, int j, typename RealBase<T>::type u );

    //------------------------------------------------------------------------//
    // Routines specific to [* ,* ] distribution                              //
    //------------------------------------------------------------------------//

    //
    // Non-collective routines
    //

    // (empty)

    //
    // Collective routines
    //
    
    // The following are all no-ops that exist to allow for more flexible 
    // templating over distribution parameters.
    template<typename S,Distribution U,Distribution V>
    void AlignWith( const DistMatrix<S,U,V>& A ) { }
    template<typename S,Distribution U,Distribution V>
    void AlignColsWith( const DistMatrix<S,U,V>& A ) { }
    template<typename S,Distribution U,Distribution V>
    void AlignRowsWith( const DistMatrix<S,U,V>& A ) { }

    // (Immutable) view of a distributed matrix
    void View( DistMatrix<T,STAR,STAR>& A );
    void LockedView( const DistMatrix<T,STAR,STAR>& A );

    // (Immutable) view of a distributed matrix's buffer
    // Create a 0 x 0 distributed matrix using the default grid
    void View
    ( int height, int width,
      T* buffer, int ldim, const elemental::Grid& grid );
    void LockedView
    ( int height, int width, 
      const T* buffer, int ldim, const elemental::Grid& grid );

    // (Immutable) view of a portion of a distributed matrix
    void View
    ( DistMatrix<T,STAR,STAR>& A, int i, int j, int height, int width );
    void LockedView
    ( const DistMatrix<T,STAR,STAR>& A, int i, int j, int height, int width );

    // (Immutable) view of two horizontally contiguous partitions of a 
    // distributed matrix
    void View1x2
    ( DistMatrix<T,STAR,STAR>& AL, DistMatrix<T,STAR,STAR>& AR );
    void LockedView1x2
    ( const DistMatrix<T,STAR,STAR>& AL, const DistMatrix<T,STAR,STAR>& AR );

    // (Immutable) view of two vertically contiguous partitions of a 
    // distributed matrix
    void View2x1
    ( DistMatrix<T,STAR,STAR>& AT,
      DistMatrix<T,STAR,STAR>& AB );
    void LockedView2x1
    ( const DistMatrix<T,STAR,STAR>& AT,
      const DistMatrix<T,STAR,STAR>& AB );

    // (Immutable) view of a contiguous 2x2 set of partitions of a 
    // distributed matrix
    void View2x2
    ( DistMatrix<T,STAR,STAR>& ATL, DistMatrix<T,STAR,STAR>& ATR,
      DistMatrix<T,STAR,STAR>& ABL, DistMatrix<T,STAR,STAR>& ABR );
    void LockedView2x2
    ( const DistMatrix<T,STAR,STAR>& ATL, const DistMatrix<T,STAR,STAR>& ATR,
      const DistMatrix<T,STAR,STAR>& ABL, const DistMatrix<T,STAR,STAR>& ABR );

    void SumOverCol();
    void SumOverRow();
    void SumOverGrid(); 

private:
    virtual void PrintBase( std::ostream& os, const std::string msg="" ) const;

    // The remainder of the class definition makes use of an idiom that allows
    // for implementing certain routines for (potentially) only complex
    // datatypes.

    template<typename Z>
    struct SetToRandomHermitianHelper
    {
        static void Func( DistMatrix<Z,STAR,STAR>& parent );
    };
    template<typename Z>
    struct SetToRandomHermitianHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,STAR>& parent );
    };
    template<typename Z> friend struct SetToRandomHermitianHelper;

    template<typename Z>
    struct SetToRandomHPDHelper
    {
        static void Func( DistMatrix<Z,STAR,STAR>& parent );
    };
    template<typename Z>
    struct SetToRandomHPDHelper<std::complex<Z> >
    {
        static void Func( DistMatrix<std::complex<Z>,STAR,STAR>& parent );
    };
    template<typename Z> friend struct SetToRandomHPDHelper;

    template<typename Z>
    struct GetRealHelper
    {
        static Z Func( const DistMatrix<Z,STAR,STAR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetRealHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,STAR>& parent, int i, int j );
    };
    template<typename Z> friend struct GetRealHelper;

    template<typename Z>
    struct GetImagHelper
    {
        static Z Func( const DistMatrix<Z,STAR,STAR>& parent, int i, int j );
    };
    template<typename Z>
    struct GetImagHelper<std::complex<Z> >
    {
        static Z Func
        ( const DistMatrix<std::complex<Z>,STAR,STAR>& parent, int i, int j ); 
    };
    template<typename Z> friend struct GetImagHelper;

    template<typename Z>
    struct SetRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,STAR>& parent, 
          int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetRealHelper;

    template<typename Z>
    struct SetImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct SetImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,STAR>& parent, 
          int i, int j, Z alpha );
    };
    template<typename Z> friend struct SetImagHelper;

    template<typename Z>
    struct UpdateRealHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateRealHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,STAR>& parent, 
          int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateRealHelper;

    template<typename Z>
    struct UpdateImagHelper
    {
        static void Func
        ( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha );
    };
    template<typename Z>
    struct UpdateImagHelper<std::complex<Z> >
    {
        static void Func
        ( DistMatrix<std::complex<Z>,STAR,STAR>& parent, 
          int i, int j, Z alpha );
    };
    template<typename Z> friend struct UpdateImagHelper;
};

} // namespace elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#include "./star_star_main.hpp"
#include "./star_star_helpers.hpp"

namespace elemental {

template<typename T>
inline
DistMatrix<T,STAR,STAR>::DistMatrix( const elemental::Grid& g )
: AbstractDistMatrix<T>
  (0,0,false,false,0,0,0,0,0,0,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,STAR>::DistMatrix
( int height, int width, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,0,0,height,width,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,STAR>::DistMatrix
( int height, int width, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,false,false,0,0,0,0,height,width,ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,STAR>::DistMatrix
( int height, int width, const T* buffer, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,0,0,0,height,width,buffer,ldim,g)
{ }

template<typename T>
inline
DistMatrix<T,STAR,STAR>::DistMatrix
( int height, int width, T* buffer, int ldim, const elemental::Grid& g )
: AbstractDistMatrix<T>
  (height,width,0,0,0,0,height,width,buffer,ldim,g)
{ }

template<typename T>
template<Distribution U,Distribution V>
inline
DistMatrix<T,STAR,STAR>::DistMatrix( const DistMatrix<T,U,V>& A )
: AbstractDistMatrix<T>(0,0,false,false,0,0,0,0,0,0,A.Grid())
{
#ifndef RELEASE
    PushCallStack("DistMatrix[* ,* ]::DistMatrix");
#endif
    if( STAR != U || STAR != V || 
        reinterpret_cast<const DistMatrix<T,STAR,STAR>*>(&A) != this )    
        *this = A;
    else
        throw std::logic_error("Tried to construct [* ,* ] with itself");
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline
DistMatrix<T,STAR,STAR>::~DistMatrix()
{ }

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetGrid( const elemental::Grid& grid )
{
    this->Empty();
    this->grid_ = &grid;
}

//
// The remainder of this file is for implementing the helpers
//

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetToRandomHermitian()
{ SetToRandomHermitianHelper<T>::Func( *this ); }

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetToRandomHPD()
{ SetToRandomHPDHelper<T>::Func( *this ); }

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,STAR>::GetReal( int i, int j ) const
{ return GetRealHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR>::GetRealHelper<Z>::Func
( const DistMatrix<Z,STAR,STAR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetRealHelper");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline typename RealBase<T>::type
DistMatrix<T,STAR,STAR>::GetImag( int i, int j ) const
{ return GetImagHelper<T>::Func( *this, i, j ); }

template<typename T>
template<typename Z>
inline Z
DistMatrix<T,STAR,STAR>::GetImagHelper<Z>::Func
( const DistMatrix<Z,STAR,STAR>& parent, int i, int j )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::GetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetReal
( int i, int j, typename RealBase<T>::type alpha )
{ SetRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetRealHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::SetImag
( int i, int j, typename RealBase<T>::type alpha )
{ SetImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::SetImagHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::SetImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::UpdateReal
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateRealHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::UpdateRealHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateReal");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

template<typename T>
inline void
DistMatrix<T,STAR,STAR>::UpdateImag
( int i, int j, typename RealBase<T>::type alpha )
{ UpdateImagHelper<T>::Func( *this, i, j, alpha ); }

template<typename T>
template<typename Z>
inline void
DistMatrix<T,STAR,STAR>::UpdateImagHelper<Z>::Func
( DistMatrix<Z,STAR,STAR>& parent, int i, int j, Z alpha )
{
#ifndef RELEASE
    PushCallStack("[* ,* ]::UpdateImag");
#endif
    throw std::logic_error("Called complex-only routine with real datatype");
}

} // elemental

#endif /* ELEMENTAL_DIST_MATRIX_STAR_STAR_HPP */

