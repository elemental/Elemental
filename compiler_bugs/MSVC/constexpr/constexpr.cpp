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

template<Dist U> constexpr Dist Collect() { return (U == CIRC) ? CIRC : STAR; }

template<typename T,Dist U,Dist V>
class DistMatrix { };

template<typename T,Dist U,Dist V>
void AllGather
(const DistMatrix<T,U,V>& A, DistMatrix<T,Collect<U>(),Collect<V>()>& B)
{ }

#ifdef USE_CONSTEXPR
template void AllGather(const DistMatrix<int,MC,MR>& A, DistMatrix<int,Collect<MC>(),Collect<MR>()>& B);
template void AllGather(const DistMatrix<double,CIRC,CIRC>& A, DistMatrix<double,Collect<CIRC>(),Collect<CIRC>()>& B);
#else
template void AllGather(const DistMatrix<int,MC,MR>& A, DistMatrix<int,STAR,STAR>& B);
template void AllGather(const DistMatrix<double,CIRC,CIRC>& A, DistMatrix<double,CIRC,CIRC>& B);
#endif

int main(int argc, char* argv[])
{
	return 0;
}
