/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SORTEIG_HPP
#define LAPACK_SORTEIG_HPP

namespace elem {

namespace internal {

template<typename R>
struct IndexValuePair {
    int index;
    R value;    

    static bool Compare
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return a.value < b.value; }
};

template<typename R>
struct IndexValuePair<Complex<R> > {
    int index;
    Complex<R> value;    

    static bool Compare
    ( const IndexValuePair<R>& a, const IndexValuePair<R>& b )
    { return Abs(a.value) < Abs(b.value); }
};

} // namespace internal

template<typename R>
inline void
SortEig( Matrix<R>& w )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    R* wBuffer = w.Buffer();
    const int k = w.Height();
    std::sort( &wBuffer[0], &wBuffer[k] );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SortEig( DistMatrix<R,VR,STAR>& w )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    const int k = w.Height();

    // Gather a full copy of w on each process and locally sort
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );
    R* wBuffer = w_STAR_STAR.LocalBuffer();
    std::sort( &wBuffer[0], &wBuffer[k] );

    // Refill the distributed w with the sorted values
    w = w_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SortEig( Matrix<R>& w, Matrix<R>& Z )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    const int n = Z.Height();
    const int k = Z.Width();

    // Initialize the pairs of indices and eigenvalues
    std::vector<internal::IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w.Get(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    std::sort
    ( pairs.begin(), pairs.end(), internal::IndexValuePair<R>::Compare );

    // Reorder the eigenvectors and eigenvalues using the new ordering
    Matrix<R> ZPerm( n, k );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), n );
        w.Set(j,0,pairs[j].value);
    }
    Z = ZPerm;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<R>& Z )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    const int n = Z.Height();
    const int k = Z.Width();
    const Grid& g = Z.Grid();

    DistMatrix<R,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );

    // Initialize the pairs of indices and eigenvalues
    std::vector<internal::IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w_STAR_STAR.GetLocal(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    std::sort
    ( pairs.begin(), pairs.end(), internal::IndexValuePair<R>::Compare );

    // Locally reorder the eigenvectors and eigenvalues using the new ordering
    const int nLocal = Z_VC_STAR.LocalHeight();
    DistMatrix<R,VC,STAR> ZPerm_VC_STAR( n, k, g );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy
        ( ZPerm_VC_STAR.LocalBuffer(0,j), 
          Z_VC_STAR.LockedLocalBuffer(0,source), nLocal );
        w_STAR_STAR.SetLocal(j,0,pairs[j].value);
    }
    Z_VC_STAR.Empty();

    Z = ZPerm_VC_STAR;
    w = w_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
SortEig( Matrix<R>& w, Matrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    const int n = Z.Height();
    const int k = Z.Width();

    // Initialize the pairs of indices and eigenvalues
    std::vector<internal::IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w.Get(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    std::sort
    ( pairs.begin(), pairs.end(), internal::IndexValuePair<R>::Compare );

    // Reorder the eigenvectors and eigenvalues using the new ordering
    Matrix<Complex<R> > ZPerm( n, k );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy( ZPerm.Buffer(0,j), Z.LockedBuffer(0,source), n );
        w.Set(j,0,pairs[j].value);
    }
    Z = ZPerm;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
SortEig( DistMatrix<R,VR,STAR>& w, DistMatrix<Complex<R> >& Z )
{
#ifndef RELEASE
    PushCallStack("SortEig");
#endif
    const int n = Z.Height();
    const int k = Z.Width();
    const Grid& g = Z.Grid();

    DistMatrix<Complex<R>,VC,STAR> Z_VC_STAR( Z );
    DistMatrix<R,STAR,STAR> w_STAR_STAR( w );

    // Initialize the pairs of indices and eigenvalues
    std::vector<internal::IndexValuePair<R> > pairs( k );
    for( int i=0; i<k; ++i )
    {
        pairs[i].index = i;
        pairs[i].value = w_STAR_STAR.GetLocal(i,0);
    }

    // Sort the eigenvalues and simultaneously form the permutation
    std::sort
    ( pairs.begin(), pairs.end(), internal::IndexValuePair<R>::Compare );

    // Locally reorder the eigenvectors and eigenvalues using the new ordering
    const int mLocal = Z_VC_STAR.LocalHeight();
    DistMatrix<Complex<R>,VC,STAR> ZPerm_VC_STAR( n, k, g );
    for( int j=0; j<k; ++j )
    {
        const int source = pairs[j].index;
        MemCopy
        ( ZPerm_VC_STAR.LocalBuffer(0,j), 
          Z_VC_STAR.LockedLocalBuffer(0,source), mLocal );
        w_STAR_STAR.SetLocal(j,0,pairs[j].value);
    }
    Z_VC_STAR.Empty();

    Z = ZPerm_VC_STAR;
    w = w_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_SORTEIG_HPP
