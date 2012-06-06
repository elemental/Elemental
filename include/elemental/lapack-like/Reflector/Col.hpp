/*
   Copyright (c) 1992-2008 The University of Tennessee.
   All rights reserved.

   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is partially based upon the LAPACK 
   routines dlarfg.f and zlarfg.f.

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

//
// Follows the LAPACK convention of defining tau such that
//
//   H = I - tau [1; v] [1, v'],
//
// but adjoint(H) [chi; x] = [beta; 0]. 
//
// Note that the adjoint of H is applied. In the case of real data,
// H' = H, so there is no complication.
//
// On exit, chi is overwritten with beta, and x is overwritten with v.
//
// The major difference from LAPACK is in the treatment of the special case 
// of x=0, where LAPACK would put H := I, which is not a valid Householder 
// reflector. We instead follow the FLAME convention of defining H such that 
//    adjoint(H) [chi; 0] = [-chi; 0],
// which is accomplished by setting tau=2, and v=0.
//

template<typename R>
inline R
ColReflector( DistMatrix<R>& chi, DistMatrix<R>& x )
{
#ifndef RELEASE
    PushCallStack("internal::ColReflector");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Width() != 1 )
        throw std::logic_error("x must be a column vector");
    if( chi.Grid().Col() != chi.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
    if( x.Grid().Col() != x.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
#endif
    const Grid& grid = x.Grid();
    mpi::Comm colComm = grid.ColComm();
    const int gridHeight = grid.Height();
    const int gridRow = grid.Row();
    const int colAlignment = chi.ColAlignment();

    std::vector<R> localNorms(gridHeight);
    R localNorm = Nrm2( x.LockedLocalMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
    R norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );

    if( norm == 0 )
    {
        if( gridRow == colAlignment )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (R)2;
    }

    R alpha;
    if( gridRow == colAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, colAlignment, colComm );

    R beta;
    if( alpha <= 0 )
        beta = lapack::SafeNorm( alpha, norm );
    else
        beta = -lapack::SafeNorm( alpha, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count = 0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        localNorm = Nrm2( x.LockedLocalMatrix() );
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
        norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );
        if( alpha <= 0 )
            beta = lapack::SafeNorm( alpha, norm );
        else
            beta = -lapack::SafeNorm( alpha, norm );
    }

    R tau = (beta-alpha)/beta;
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridRow == colAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

template<typename R> 
inline Complex<R>
ColReflector
( DistMatrix<Complex<R> >& chi, DistMatrix<Complex<R> >& x )
{
#ifndef RELEASE
    PushCallStack("internal::ColReflector");
    if( chi.Grid() != x.Grid() )
        throw std::logic_error
        ("chi and x must be distributed over the same grid");
    if( chi.Height() != 1 || chi.Width() != 1 )
        throw std::logic_error("chi must be a scalar");
    if( x.Width() != 1 )
        throw std::logic_error("x must be a column vector");
    if( chi.Grid().Col() != chi.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
    if( x.Grid().Col() != x.RowAlignment() )
        throw std::logic_error("Reflecting with incorrect column of processes");
#endif
    typedef Complex<R> C;
    const Grid& grid = x.Grid();
    mpi::Comm colComm = grid.ColComm();
    const int gridHeight = grid.Height();
    const int gridRow = grid.Row();
    const int colAlignment = chi.ColAlignment();

    std::vector<R> localNorms(gridHeight);
    R localNorm = Nrm2( x.LockedLocalMatrix() ); 
    mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
    R norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );

    C alpha;
    if( gridRow == colAlignment )
        alpha = chi.GetLocal(0,0);
    mpi::Broadcast( &alpha, 1, colAlignment, colComm );

    if( norm == (R)0 && alpha.imag == (R)0 )
    {
        if( gridRow == colAlignment )
            chi.SetLocal(0,0,-chi.GetLocal(0,0));
#ifndef RELEASE
        PopCallStack();
#endif
        return (C)2;
    }

    R beta;
    if( alpha.real <= 0 )
        beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
    else
        beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );

    const R one = 1;
    const R safeMin = lapack::MachineSafeMin<R>();
    const R epsilon = lapack::MachineEpsilon<R>();
    const R safeInv = safeMin/epsilon;
    int count = 0;
    if( Abs(beta) < safeInv )
    {
        R invOfSafeInv = one/safeInv;
        do
        {
            ++count;
            Scale( invOfSafeInv, x );
            alpha *= invOfSafeInv;
            beta *= invOfSafeInv;
        } while( Abs(beta) < safeInv );

        localNorm = Nrm2( x.LockedLocalMatrix() );
        mpi::AllGather( &localNorm, 1, &localNorms[0], 1, colComm );
        norm = blas::Nrm2( gridHeight, &localNorms[0], 1 );
        if( alpha.real <= 0 )
            beta = lapack::SafeNorm( alpha.real, alpha.imag, norm );
        else
            beta = -lapack::SafeNorm( alpha.real, alpha.imag, norm );
    }

    C tau = C( (beta-alpha.real)/beta, -alpha.imag/beta );
    Scale( one/(alpha-beta), x );

    for( int j=0; j<count; ++j )
        beta *= safeInv;
    if( gridRow == colAlignment )
        chi.SetLocal(0,0,beta);
        
#ifndef RELEASE
    PopCallStack();
#endif
    return tau;
}

} // namespace internal
} // namespace elem
