/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

Int NaturalBisect
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph,
        Int& nxLeft,
        Int& nyLeft,
        Int& nzLeft,
        Graph& leftChild,
        Int& nxRight,
        Int& nyRight,
        Int& nzRight,
        Graph& rightChild,
        vector<Int>& perm )
{
    EL_DEBUG_CSE
    const Int numSources = graph.NumSources();
    if( numSources == 0 )
        LogicError("There is no reason to bisect an empty sequential graph");

    Int leftChildSize, rightChildSize, sepSize;
    perm.resize( numSources );
    if( nx >= ny && nx >= nz )
    {
        nxLeft = (nx-1)/2;
        nyLeft = ny;
        nzLeft = nz;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx-1-nxLeft;
        nyRight = ny;
        nzRight = nz;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = ny*nz;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nxLeft; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=nxLeft+1; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                perm[nxLeft+y*nx+z*nx*ny] = off++;
    }
    else if( ny >= nx && ny >= nz )
    {
        nxLeft = nx;
        nyLeft = (ny-1)/2;
        nzLeft = nz;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx;
        nyRight = ny-1-nyLeft;
        nzRight = nz;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = nx*nz;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nz; ++z )
            for( Int y=0; y<nyLeft; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int y=nyLeft+1; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int z=0; z<nz; ++z )
            for( Int x=0; x<nx; ++x )
                perm[x+nyLeft*nx+z*nx*ny] = off++;
    }
    else
    {
        nxLeft = nx;
        nyLeft = ny;
        nzLeft = (nz-1)/2;
        leftChildSize = nxLeft*nyLeft*nzLeft;

        nxRight = nx;
        nyRight = ny;
        nzRight = nz-1-nzLeft;
        rightChildSize = nxRight*nyRight*nzRight;

        sepSize = nx*ny;

        // Fill the left side
        Int off=0;
        for( Int z=0; z<nzLeft; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the right side
        off = leftChildSize;
        for( Int z=nzLeft+1; z<nz; ++z )
            for( Int y=0; y<ny; ++y )
                for( Int x=0; x<nx; ++x )
                    perm[x+y*nx+z*nx*ny] = off++;

        // Fill the separator
        off=leftChildSize+rightChildSize;
        for( Int y=0; y<ny; ++y )
            for( Int x=0; x<nx; ++x )
                perm[x+y*nx+nzLeft*nx*ny] = off++;
    }
    EL_DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildrenFromPerm
    ( graph, perm, leftChildSize, leftChild, rightChildSize, rightChild );
    return sepSize;
}

Int NaturalBisect
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph,
        Int& nxChild,
        Int& nyChild,
        Int& nzChild,
        unique_ptr<Grid>& childGrid,
        DistGraph& child,
        DistMap& perm,
        bool& onLeft )
{
    EL_DEBUG_CSE
    const Int numSources = graph.NumSources();
    const Int firstLocalSource = graph.FirstLocalSource();
    const Int numLocalSources = graph.NumLocalSources();
    const Grid& grid = graph.Grid();
    const int commSize = grid.Size();
    if( commSize == 1 )
        LogicError
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    Int leftChildSize, rightChildSize, sepSize;
    Int nxLeft, nyLeft, nzLeft, nxRight, nyRight, nzRight;
    perm.SetGrid( grid );
    perm.Resize( numSources );
    if( nx != 0 && ny != 0 && nz != 0 )
    {
        if( nx >= ny && nx >= nz )
        {
            nxLeft = (nx-1)/2;
            nyLeft = ny;
            nzLeft = nz;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx-1-nxLeft;
            nyRight = ny;
            nzRight = nz;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = ny*nz;

            const Int rightOff=leftChildSize,
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( x < nxLeft )
                {
                    const Int xLeft = x;
                    const Int leftInd = xLeft + y*nxLeft + z*nxLeft*ny;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( x > nxLeft )
                {
                    const Int xRight = x-(nxLeft+1);
                    const Int rightInd = xRight + y*nxRight + z*nxRight*ny;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = y + z*ny;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
        else if( ny >= nx && ny >= nz )
        {
            nxLeft = nx;
            nyLeft = (ny-1)/2;
            nzLeft = nz;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx;
            nyRight = ny-1-nyLeft;
            nzRight = nz;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = nx*nz;

            const Int rightOff=leftChildSize,
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( y < nyLeft )
                {
                    const Int yLeft = y;
                    const Int leftInd = x + yLeft*nx + z*nx*nyLeft;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( y > nyLeft )
                {
                    const Int yRight = y - (nyLeft+1);
                    const Int rightInd = x + yRight*nx + z*nx*nyRight;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = x + z*nx;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
        else
        {
            nxLeft = nx;
            nyLeft = ny;
            nzLeft = (nz-1)/2;
            leftChildSize = nxLeft*nyLeft*nzLeft;

            nxRight = nx;
            nyRight = ny;
            nzRight = nz-1-nzLeft;
            rightChildSize = nxRight*nyRight*nzRight;

            sepSize = nx*ny;

            const Int rightOff=leftChildSize,
                      sepOff=leftChildSize+rightChildSize;
            for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
            {
                const Int i = iLocal + firstLocalSource;
                const Int x = i % nx;
                const Int y = (i/nx) % ny;
                const Int z = i/(nx*ny);
                if( z < nzLeft )
                {
                    const Int zLeft = z;
                    const Int leftInd = x + y*nx + zLeft*nx*ny;
                    perm.SetLocal( iLocal, leftInd );
                }
                else if( z > nzLeft )
                {
                    const Int zRight = z - (nzLeft+1);
                    const Int rightInd = x + y*nx + zRight*nx*ny;
                    perm.SetLocal( iLocal, rightOff+rightInd );
                }
                else
                {
                    const Int sepInd = x + y*nx;
                    perm.SetLocal( iLocal, sepOff+sepInd );
                }
            }
        }
    }
    else
    {
        leftChildSize = rightChildSize = sepSize = 0;
        nxLeft = nx;
        nyLeft = ny;
        nzLeft = nz;
        nxRight = nx;
        nyRight = ny;
        nzRight = nz;
    }
    EL_DEBUG_ONLY(EnsurePermutation( perm ))

    BuildChildFromPerm
    ( graph, perm, leftChildSize, rightChildSize, onLeft, childGrid, child );

    if( onLeft )
    {
        nxChild = nxLeft;
        nyChild = nyLeft;
        nzChild = nzLeft;
    }
    else
    {
        nxChild = nxRight;
        nyChild = nyRight;
        nzChild = nzRight;
    }
    return sepSize;
}

} // namespace El
