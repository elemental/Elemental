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

inline void
elemental::advanced::internal::ComposePivots
( const DistMatrix<int,STAR,STAR>& p, 
  std::vector<int>& image, std::vector<int>& preimage, int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ComposePivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( pivotOffset < 0 )
        throw std::logic_error("pivotOffset must be non-negative");
#endif
    const int b = p.Height();
    image.resize( b );
    preimage.resize( b );

    // Construct the image of {0,...,b-1} under the permutation
    for( int i=0; i<b; ++i )
    {
        int row = i;
        for( int j=0; j<min(b,row+1); ++j )
        {
            if( p.GetLocalEntry(j,0)-pivotOffset == row )
                row = j;
            else if( j == row )
                row = p.GetLocalEntry(j,0)-pivotOffset;
        }
        image[i] = row;
    }

    // Construct the preimage of {0,...,b-1} under the permutation
    for( int i=0; i<b; ++i )
    {
        int row = p.GetLocalEntry(i,0)-pivotOffset;
        for( int j=i-1; j>=0; --j )
        {
            if( p.GetLocalEntry(j,0)-pivotOffset == row )
                row = j;
            else if( j == row )
                row = p.GetLocalEntry(j,0)-pivotOffset;
        }
        preimage[i] = row;
    }

#ifndef RELEASE
    PopCallStack();
#endif
}

