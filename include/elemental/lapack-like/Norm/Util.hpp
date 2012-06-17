/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename F,Distribution U,Distribution V>
mpi::Comm 
NormComm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::NormComm");
#endif
    const Grid& grid = A.Grid();
    mpi::Comm comm;
    if( U == MC && V == MR )
        comm = grid.VCComm();
    else if( U == MC && V == STAR )
        comm = grid.MCComm();
    else if( U == MD && V == STAR )
        comm = grid.VCComm();
    else if( U == MR && V == MC )
        comm = grid.VRComm();
    else if( U == MR && V == STAR )
        comm = grid.MRComm();
    else if( U == STAR && V == MC )
        comm = grid.MCComm();
    else if( U == STAR && V == MD )
        comm = grid.VCComm();
    else if( U == STAR && V == MR )
        comm = grid.MRComm();
    else if( U == STAR && V == STAR )
        comm = mpi::COMM_SELF;
    else if( U == STAR && V == VC )
        comm = grid.VCComm();
    else if( U == STAR && V == VR )
        comm = grid.VRComm();
    else if( U == VC && V == STAR )
        comm = grid.VCComm();
    else if( U == VR && V == STAR )
        comm = grid.VRComm();
    else
        throw std::logic_error("Invalid distribution");
#ifndef RELEASE
    PopCallStack();
#endif
    return comm;
}

template<typename F,Distribution U,Distribution V>
mpi::Comm 
NormColComm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::NormComm");
#endif
    const Grid& grid = A.Grid();
    mpi::Comm comm;
    switch( U )
    {
    case MC: comm = grid.MCComm(); break;
    case MD: comm = grid.VCComm(); break;
    case MR: comm = grid.MRComm(); break;
    case VC: comm = grid.VCComm(); break;
    case VR: comm = grid.VRComm(); break;
    case STAR: comm = mpi::COMM_SELF; break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return comm;
}

template<typename F,Distribution U,Distribution V>
mpi::Comm 
NormRowComm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("internal::NormComm");
#endif
    const Grid& grid = A.Grid();
    mpi::Comm comm;
    switch( V )
    {
    case MC: comm = grid.MCComm(); break;
    case MD: comm = grid.VCComm(); break;
    case MR: comm = grid.MRComm(); break;
    case VC: comm = grid.VCComm(); break;
    case VR: comm = grid.VRComm(); break;
    case STAR: comm = mpi::COMM_SELF; break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return comm;
}

} // namespace internal
} // namespace elem
