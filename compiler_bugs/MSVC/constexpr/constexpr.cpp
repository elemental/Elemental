/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <iostream>

/* The following is meant to test for what appears to be a bug in the handling
   of constexpr's in template instantiations in Visual Studio 2015 (pre-release)
*/

enum Dist
{
  MC,
  MD,
  MR,
  VC,
  VR,
  STAR,
  CIRC
};

template<Dist U> constexpr Dist Collect()       { return STAR; }
template<>       constexpr Dist Collect<CIRC>() { return CIRC; }

template<typename T,Dist U,Dist V>
class DistMatrix { };

template<typename T,Dist U,Dist V>
void AllGather
( const DistMatrix<T,        U,           V   >& A, 
        DistMatrix<T,Collect<U>(),Collect<V>()>& B ) 
{ 
    std::cout << "U=" << U << ", V=" << V << std::endl; 
}

int main( int argc, char* argv[] )
{
    {
        DistMatrix<double,MC,MR> A;
        DistMatrix<double,STAR,STAR> B;
        AllGather( A, B );
    }
    {
        DistMatrix<int,CIRC,CIRC> A;
        DistMatrix<int,CIRC,CIRC> B;
        AllGather( A, B );
    }
    return 0;
}
