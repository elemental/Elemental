/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FRONTUTIL_HPP
#define EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FRONTUTIL_HPP

namespace El {
namespace ldl {

template<typename F>
void FormDiagonalBlocks
( const DistMatrix<F,VC,STAR>& L, DistMatrix<F,STAR,STAR>& D, bool conjugate )
{
    const Int height = L.Width();
    auto LT = L( IR(0,height), IR(0,height) );

    const Int blocksize = Blocksize();
    const int colStride = LT.ColStride();
    const Int localHeight = LT.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int portionSize = maxLocalHeight*blocksize;

    vector<F> sendBuf( portionSize );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = L.GlobalRow(iLoc);
        const Int block = i / blocksize;
        const Int jStart = block*blocksize;
        const Int b = Min(height-jStart,blocksize);
        for( Int jOff=0; jOff<b; ++jOff )
        {
            const F val = L.GetLocal(iLoc,jStart+jOff);
            sendBuf[iLoc*blocksize+jOff] = ( conjugate ? Conj(val) : val );
        }
    }

    vector<F> recvBuf( portionSize*colStride );
    mpi::AllGather
    ( sendBuf.data(), portionSize, 
      recvBuf.data(), portionSize, L.ColComm() );
    SwapClear( sendBuf );
    
    D.Resize( blocksize, height );
    for( Int q=0; q<colStride; ++q )
    {
        const F* pkg = &recvBuf[q*portionSize];
        const Int qLocalHeight = Length(height,q,colStride);
        for( Int iLoc=0; iLoc<qLocalHeight; ++iLoc )
        {
            const Int i = q + iLoc*colStride;
            for( Int jOff=0; jOff<blocksize; ++jOff )
                D.SetLocal( jOff, i, pkg[jOff+iLoc*blocksize] );
        }
    }
}

template<typename F>
void AccumulateRHS( const DistMatrix<F,VC,STAR>& X, DistMatrix<F,STAR,STAR>& Z )
{
    const Int height = X.Height();
    const Int width = X.Width();
    Z.Empty();
    Zeros( Z, height, width );

    const Int localHeight = X.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = X.GlobalRow(iLoc);
        for( Int j=0; j<width; ++j )
            Z.SetLocal( i, j, X.GetLocal(iLoc,j) );
    }
    El::AllReduce( Z, X.ColComm() );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_NUMERIC_LOWERSOLVE_FRONTUTIL_HPP
