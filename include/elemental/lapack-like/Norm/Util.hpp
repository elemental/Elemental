/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F,Distribution U,Distribution V>
inline mpi::Comm 
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
inline mpi::Comm 
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
inline mpi::Comm 
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
