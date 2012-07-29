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

//----------------------------------------------------------------------------//
// Local LAPACK-like routines                                                 //
//----------------------------------------------------------------------------//

template<typename F>
void LocalCholesky( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalHegst
( LeftOrRight side, UpperOrLower uplo, 
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B );

template<typename F>
void LocalLDL
( Orientation orientation, 
  DistMatrix<F,STAR,STAR>& A, DistMatrix<F,STAR,STAR>& d );

template<typename F>
void LocalLU( DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

template<typename F>
void LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A );

//----------------------------------------------------------------------------//
// Cholesky helpers                                                           //
//----------------------------------------------------------------------------//

template<typename F>
void CholeskyLVar2( DistMatrix<F>& A );

template<typename F>
void CholeskyLVar2Naive( DistMatrix<F>& A );

template<typename F>
void CholeskyLVar3( DistMatrix<F>& A );

template<typename F>
void CholeskyLVar3Naive( DistMatrix<F>& A );

template<typename F>
void CholeskyLVar3Square( DistMatrix<F>& A );

template<typename F>
void CholeskyUVar2( DistMatrix<F>& A );

template<typename F>
void CholeskyUVar2Naive( DistMatrix<F>& A );
 
template<typename F>
void CholeskyUVar3( DistMatrix<F>& A );

template<typename F>
void CholeskyUVar3Naive( DistMatrix<F>& A );

template<typename F>
void CholeskyUVar3Square( DistMatrix<F>& A );
            
//----------------------------------------------------------------------------//
// GaussElim                                                                  //
//----------------------------------------------------------------------------//
            
template<typename F>
void ReduceToRowEchelon( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void ReduceToRowEchelon( DistMatrix<F>& A, DistMatrix<F>& B );

//----------------------------------------------------------------------------//
// Hegst                                                                      //
//----------------------------------------------------------------------------//

template<typename F>
void HegstLLVar1( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstLLVar2( DistMatrix<F>& A, const DistMatrix<F>& L );

// HegstLLVar3 would redundantly compute too much data

template<typename F>
void HegstLLVar4( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstLLVar5( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstLUVar1( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstLUVar2( DistMatrix<F>& A, const DistMatrix<F>& U );

// HegstLUVar3 would redundantly compute too much data

template<typename F>
void HegstLUVar4( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstLUVar5( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstRLVar1( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstRLVar2( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstRLVar3( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstRLVar4( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstRLVar5( DistMatrix<F>& A, const DistMatrix<F>& L );

template<typename F>
void HegstRUVar1( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstRUVar2( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstRUVar3( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstRUVar4( DistMatrix<F>& A, const DistMatrix<F>& U );

template<typename F>
void HegstRUVar5( DistMatrix<F>& A, const DistMatrix<F>& U );

//----------------------------------------------------------------------------//
// LDL                                                                        //
//----------------------------------------------------------------------------//

template<typename F>
void LDLVar3
( Orientation orientation, Matrix<F>& A, Matrix<F>& d );

template<typename F>
void LDLVar3
( Orientation orientation, DistMatrix<F>& A, DistMatrix<F,MC,STAR>& d );

//----------------------------------------------------------------------------//
// LU                                                                         //
//----------------------------------------------------------------------------//

void ComposePanelPivots
( const Matrix<int>& p,
        int pivotOffset,
        std::vector<int>& image,
        std::vector<int>& preimage );

void ComposePanelPivots
( const DistMatrix<int,STAR,STAR>& p,
        int pivotOffset,
        std::vector<int>& image,
        std::vector<int>& preimage );

template<typename F> void CreatePivotOp();
template<> void CreatePivotOp<float>();
template<> void CreatePivotOp<double>();
template<> void CreatePivotOp<scomplex>();
template<> void CreatePivotOp<dcomplex>();

template<typename T> void DestroyPivotOp();
template<> void DestroyPivotOp<float>();
template<> void DestroyPivotOp<double>();
template<> void DestroyPivotOp<scomplex>();
template<> void DestroyPivotOp<dcomplex>();

template<typename F>
void PanelLU( Matrix<F>& A, Matrix<int>& p, int pivotOffset=0 );
template<typename F>
void PanelLU
( DistMatrix<F,STAR,STAR>& A, 
  DistMatrix<F,MC,  STAR>& B, 
  DistMatrix<int,STAR,STAR>& p, 
  int pivotOffset=0 );

template<typename F> mpi::Op PivotOp();
template<> mpi::Op PivotOp<float>();
template<> mpi::Op PivotOp<double>();
template<> mpi::Op PivotOp<scomplex>();
template<> mpi::Op PivotOp<dcomplex>();
            
template<typename F>
void PivotFunc
( void* inData, void* outData, int* length, mpi::Datatype* datatype );

//----------------------------------------------------------------------------//
// LQ                                                                         //
//----------------------------------------------------------------------------//

template<typename R>
void PanelLQ( Matrix<R>& A );
template<typename R>
void PanelLQ( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );

template<typename R>
void PanelLQ( DistMatrix<R>& A );
template<typename R>
void PanelLQ
( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t );

//----------------------------------------------------------------------------//
// QR                                                                         //
//----------------------------------------------------------------------------//

template<typename R>
void PanelQR( Matrix<R>& A );
template<typename R>
void PanelQR( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t );

template<typename R>
void PanelQR( DistMatrix<R>& A );
template<typename R>
void PanelQR
( DistMatrix<Complex<R> >& A, DistMatrix<Complex<R>,MD,STAR>& t );

//----------------------------------------------------------------------------//
// Reflector                                                                  //
//----------------------------------------------------------------------------//
template<typename R>
R ColReflector( DistMatrix<R>& chi, DistMatrix<R>& x );

template<typename R>
Complex<R> ColReflector
( DistMatrix<Complex<R> >& chi, 
  DistMatrix<Complex<R> >& x );

template<typename R>
R RowReflector( DistMatrix<R>& chi, DistMatrix<R>& x );

template<typename R>
Complex<R> RowReflector
( DistMatrix<Complex<R> >& chi,
  DistMatrix<Complex<R> >& x );

//----------------------------------------------------------------------------//
// Bidiag                                                                     //
//----------------------------------------------------------------------------//

template<typename R>
void BidiagL( DistMatrix<R>& A );
template<typename R>
void BidiagU( DistMatrix<R>& A );

template<typename R>
void BidiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ );
template<typename R>
void BidiagU
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ );

template<typename R>
void PanelBidiagL
( DistMatrix<R>& A,
  DistMatrix<R>& X,
  DistMatrix<R>& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR );
template<typename R>
void PanelBidiagU
( DistMatrix<R>& A,
  DistMatrix<R>& X,
  DistMatrix<R>& Y,
  DistMatrix<R,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<R,STAR,MR  >& ARowPan_STAR_MR );

template<typename R>
void PanelBidiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,MD,  STAR>& tP,
  DistMatrix<Complex<R>,MD,  STAR>& tQ,
  DistMatrix<Complex<R> >& X,
  DistMatrix<Complex<R> >& Y,
  DistMatrix<Complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<Complex<R>,STAR,MR  >& ARowPan_STAR_MR );
template<typename R>
void PanelBidiagU
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,MD,  STAR>& tP,
  DistMatrix<Complex<R>,MD,  STAR>& tQ,
  DistMatrix<Complex<R> >& X,
  DistMatrix<Complex<R> >& Y,
  DistMatrix<Complex<R>,MC,  STAR>& AColPan_MC_STAR,
  DistMatrix<Complex<R>,STAR,MR  >& ARowPan_STAR_MR );

//----------------------------------------------------------------------------//
// HermitianTridiag                                                           //
//----------------------------------------------------------------------------//

template<typename R>
void HermitianPanelTridiagL
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagU
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagLSquare
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagUSquare
( DistMatrix<R>& A, 
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR,
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR );

template<typename R>
void HermitianPanelTridiagL
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagU
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );

template<typename R>
void HermitianPanelTridiagLSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
template<typename R>
void HermitianPanelTridiagUSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR );
 
template<typename R>
void HermitianTridiagL( DistMatrix<R>& A );
template<typename R>
void HermitianTridiagU( DistMatrix<R>& A );

template<typename R>
void HermitianTridiagLSquare( DistMatrix<R>& A );
template<typename R>
void HermitianTridiagUSquare( DistMatrix<R>& A );

template<typename R>
void HermitianTridiagL
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );
template<typename R>
void HermitianTridiagU
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );

template<typename R>
void HermitianTridiagLSquare
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );
template<typename R>
void HermitianTridiagUSquare
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,STAR,STAR>& t );

//----------------------------------------------------------------------------//
// HPD Inverse                                                                //
//----------------------------------------------------------------------------//

template<typename F>
void HPDInverseLVar2( DistMatrix<F>& A );

template<typename F>
void HPDInverseUVar2( DistMatrix<F>& A );

//----------------------------------------------------------------------------//
// Triangular Inverse                                                         //
//----------------------------------------------------------------------------//

template<typename F>
void TriangularInverseVar3
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F>& A  );

template<typename F>
void TriangularInverseLVar3
( UnitOrNonUnit diag, DistMatrix<F>& L );

template<typename F>
void TriangularInverseUVar3
( UnitOrNonUnit diag, DistMatrix<F>& U );

//----------------------------------------------------------------------------//
// LAPACK-like Utility Functions                                              //
//----------------------------------------------------------------------------//
template<typename F>
double CholeskyGFlops( int m, double seconds );

template<typename F>
double HegstGFlops( int m, double seconds );

template<typename F>
double LDLHGFlops( int m, double seconds );

template<typename F>
double LDLTGFlops( int m, double seconds );

template<typename F>
double LUGFlops( int m, double seconds );

template<typename F>
double LQGFlops( int m, int n, double seconds );

template<typename F>
double QRGFlops( int m, int n, double seconds );

template<typename F>
double HermitianTridiagGFlops( int m, double seconds );

template<typename F>
double HPDInverseGFlops( int m, double seconds );

template<typename F>
double TriangularInverseGFlops( int m, double seconds );

template<typename F>
double ApplyPackedReflectorsGFlops( int m, double seconds );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// Local LAPACK-like routines
//

template<typename F>
inline void
LocalCholesky
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::LocalCholesky");
#endif
    Cholesky( uplo, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalHegst
( LeftOrRight side, UpperOrLower uplo,
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B )
{
#ifndef RELEASE
    PushCallStack("internal::LocalHegst");
#endif
    Hegst( side, uplo, A.LocalMatrix(), B.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalLDL
( Orientation orientation, 
  DistMatrix<F,STAR,STAR>& A, DistMatrix<F,STAR,STAR>& d )
{
#ifndef RELEASE
    PushCallStack("internal::LocalLDL");
    if( d.Viewing() && (d.Height() != A.Height() || d.Width() != 1) )
        throw std::logic_error
        ("d must be a column vector of the same height as A");
#endif
    if( !d.Viewing() )
        d.ResizeTo( A.Height(), 1 );
    internal::LDLVar3
    ( orientation, A.LocalMatrix(), d.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalLU( DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("internal::LocalLU");
#endif
    LU( A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("internal::LocalHPDInverse");
#endif
    HPDInverse( uplo, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalTriangularInverse
( UpperOrLower uplo, UnitOrNonUnit diag, DistMatrix<F,STAR,STAR>& A )
{ 
#ifndef RELEASE
    PushCallStack("internal::LocalTriangularInverse");
#endif
    TriangularInverse( uplo, diag, A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// GFlop helpers
//

template<>
inline double
CholeskyGFlops<float>
( int m, double seconds )
{ return (1./3.*m*m*m)/(1.e9*seconds); }
            
template<>
inline double
CholeskyGFlops<double>
( int m, double seconds )
{ return CholeskyGFlops<float>(m,seconds); }
            
template<>
inline double
CholeskyGFlops<scomplex>
( int m, double seconds )
{ return 4.*CholeskyGFlops<float>(m,seconds); }
            
template<>
inline double
CholeskyGFlops<dcomplex>
( int m, double seconds )
{ return 4.*CholeskyGFlops<float>(m,seconds); }

template<>
inline double
HegstGFlops<float>
( int m, double seconds )
{ return (1.*m*m*m)/(1.e9*seconds); }

template<>
inline double
HegstGFlops<double>
( int m, double seconds )
{ return HegstGFlops<float>(m,seconds); }

template<>
inline double
HegstGFlops<scomplex>
( int m, double seconds )
{ return 4.*HegstGFlops<float>(m,seconds); }

template<>
inline double
HegstGFlops<dcomplex>
( int m, double seconds )
{ return 4.*HegstGFlops<float>(m,seconds); }

template<>
inline double
LUGFlops<float>
( int m, double seconds )
{ return (2./3.*m*m*m)/(1.e9*seconds); }

template<typename F>
inline double LDLHGFlops( int n, double seconds )
{ return CholeskyGFlops<F>( n, seconds ); }

template<typename F>
inline double LDLTGFlops( int n, double seconds )
{ return CholeskyGFlops<F>( n, seconds ); }

template<>
inline double
LUGFlops<double>
( int m, double seconds )
{ return LUGFlops<float>(m,seconds); }

template<>
inline double
LUGFlops<scomplex>
( int m, double seconds )
{ return 4.*LUGFlops<float>(m,seconds); }

template<>
inline double
LUGFlops<dcomplex>
( int m, double seconds )
{ return 4.*LUGFlops<float>(m,seconds); }

template<>
inline double
LQGFlops<float>
( int m, int n, double seconds )
{ return (2.*m*m*n-2./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
LQGFlops<double>
( int m, int n, double seconds )
{ return LQGFlops<float>(m,n,seconds); }

template<>
inline double
LQGFlops<scomplex>
( int m, int n, double seconds )
{ return 4.*LQGFlops<float>(m,n,seconds); }

template<>
inline double
LQGFlops<dcomplex>
( int m, int n, double seconds )
{ return 4.*LQGFlops<float>(m,n,seconds); }

template<>
inline double
QRGFlops<float>
( int m, int n, double seconds )
{ return (2.*m*n*n-2./3.*n*n*n)/(1.e9*seconds); }

template<>
inline double
QRGFlops<double>
( int m, int n, double seconds )
{ return QRGFlops<float>(m,n,seconds); }

template<>
inline double
QRGFlops<scomplex>
( int m, int n, double seconds )
{ return 4.*QRGFlops<float>(m,n,seconds); }

template<>
inline double
QRGFlops<dcomplex>
( int m, int n, double seconds )
{ return 4.*QRGFlops<float>(m,n,seconds); }

template<>
inline double
HermitianTridiagGFlops<float>
( int m, double seconds )
{ return (4./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
HermitianTridiagGFlops<double>
( int m, double seconds )
{ return HermitianTridiagGFlops<float>(m,seconds); }

template<>
inline double
HermitianTridiagGFlops<scomplex>
( int m, double seconds )
{ return 4.*HermitianTridiagGFlops<float>(m,seconds); }

template<>
inline double
HermitianTridiagGFlops<dcomplex>
( int m, double seconds )
{ return 4.*HermitianTridiagGFlops<float>(m,seconds); }

template<>
inline double
HPDInverseGFlops<float>
( int m, double seconds )
{ return (1.*m*m*m)/(1.e9*seconds); }

template<>
inline double
HPDInverseGFlops<double>
( int m, double seconds )
{ return HPDInverseGFlops<float>(m,seconds); }

template<>
inline double
HPDInverseGFlops<scomplex>
( int m, double seconds )
{ return 4.*HPDInverseGFlops<float>(m,seconds); }

template<>
inline double
HPDInverseGFlops<dcomplex>
( int m, double seconds )
{ return 4.*HPDInverseGFlops<float>(m,seconds); }

template<>
inline double
TriangularInverseGFlops<float>
( int m, double seconds )
{ return (1./3.*m*m*m)/(1.e9*seconds); }

template<>
inline double
TriangularInverseGFlops<double>
( int m, double seconds )
{ return TriangularInverseGFlops<float>(m,seconds); }

template<>
inline double
TriangularInverseGFlops<scomplex>
( int m, double seconds )
{ return 4.*TriangularInverseGFlops<float>(m,seconds); }

template<>
inline double
TriangularInverseGFlops<dcomplex>
( int m, double seconds )
{ return 4.*TriangularInverseGFlops<float>(m,seconds); }

template<>
inline double
ApplyPackedReflectorsGFlops<float>
( int m, double seconds )
{ return (2.*m*m*m)/(1.e9*seconds); }

template<>
inline double
ApplyPackedReflectorsGFlops<double>
( int m, double seconds )
{ return ApplyPackedReflectorsGFlops<float>(m,seconds); }

template<>
inline double
ApplyPackedReflectorsGFlops<scomplex>
( int m, double seconds )
{ return 4.*ApplyPackedReflectorsGFlops<float>(m,seconds); }

template<>
inline double
ApplyPackedReflectorsGFlops<dcomplex>
( int m, double seconds )
{ return 4.*ApplyPackedReflectorsGFlops<float>(m,seconds); }

} // namespace internal
} // namespace elem
