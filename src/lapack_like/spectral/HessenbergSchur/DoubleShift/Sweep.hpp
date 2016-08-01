/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_DOUBLE_SHIFT_SWEEP_HPP
#define EL_SCHUR_HESS_DOUBLE_SHIFT_SWEEP_HPP

namespace El {
namespace hess_schur {
namespace double_shift {

template<typename Real>
void PrepareShifts
( Real eta00, Real eta01,
  Real eta10, Real eta11,
  Complex<Real>& shift0,
  Complex<Real>& shift1 )
{
    DEBUG_CSE
    const Real zero(0);
    const Real scale = Abs(eta00) + Abs(eta01) + Abs(eta10) + Abs(eta11);
    if( scale == zero )
    {
        shift0 = shift1 = zero;
    }
    else
    {
        eta00 /= scale;
        eta01 /= scale;
        eta10 /= scale;
        eta11 /= scale;
        Real halfTrace = (eta00+eta11) / 2;
        Real det = (eta00-halfTrace)*(eta11-halfTrace) - eta01*eta10;
        Real absDisc = Sqrt(Abs(det));
        if( det >= zero )
        {
            shift0 = scale*Complex<Real>(halfTrace,absDisc);
            shift1 = Conj(shift0);
        }
        else
        {
            if( Abs(halfTrace+absDisc-eta11) <= Abs(halfTrace-absDisc-eta11) )
            {
                shift0 = shift1 = (halfTrace+absDisc)*scale;
            }
            else
            {
                shift0 = shift1 = (halfTrace-absDisc)*scale;
            }
        }
    }
}

template<typename Real>
Int ChooseStart
( const Matrix<Real>& H, 
  const Complex<Real>& shift0,
  const Complex<Real>& shift1,
        vector<Real>& v )
{
    DEBUG_CSE
    const Real ulp = limits::Precision<Real>();
    const Int n = H.Height();

    v.resize( 3 );
    for( Int k=n-4; k+1>=0; --k )
    {
        const Real& eta11 = H(k+1,k+1);
        const Real& eta12 = H(k+1,k+2);
        const Real& eta21 = H(k+2,k+1);
        const Real& eta22 = H(k+2,k+2);
        const Real& eta32 = H(k+3,k+2);
  
        Real scale = OneAbs(eta11-shift1) + Abs(eta21);
        Real eta21Scale = eta21 / scale;

        v[0] =
          eta21Scale*eta12 +
          (eta11-shift0.real())*((eta11-shift1.real())/scale) -
          shift0.imag()*(shift1.imag()/scale);
        v[1] = eta21Scale*(eta11+eta22-shift0.real()-shift1.real());
        v[2] = eta21Scale*eta32;

        scale = Abs(v[0]) + Abs(v[1]) + Abs(v[2]);
        v[0] /= scale;
        v[1] /= scale;
        v[2] /= scale;
        if( k == -1 )
        {
            return 0;
        }
        
        const Real& eta00 = H(k,  k);
        const Real& eta10 = H(k+1,k);
        if( Abs(eta10)*(Abs(v[1])+Abs(v[2])) <= 
            ulp*Abs(v[0])*(Abs(eta00)+Abs(eta11)+Abs(eta22)) )
        {
            return k+1;
        }
    }
    // This should never be reached but is to avoid compiler warnings
    return 0;
}

template<typename Real>
void Sweep
( Matrix<Real>& H,
  const Complex<Real>& shift0,
  const Complex<Real>& shift1,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Real zero(0), one(1);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );

    const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg ); 
    const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );

    vector<Real> v(3);
    auto subInd = IR(winBeg,winEnd);
    Int shiftStart = winBeg +
      ChooseStart( H(subInd,subInd), shift0, shift1, v );

    for( Int k=shiftStart; k<winEnd-1; ++k )
    {
        const Int numReflect = Min( 3, winEnd-k );
        if( k > shiftStart )
        {
            MemCopy( v.data(), &H(k,k-1), numReflect );
        }
        Real tau0 = lapack::Reflector( numReflect, v[0], &v[1], 1 );
        if( k > shiftStart )
        {
            H(k,  k-1) = v[0];
            H(k+1,k-1) = zero;
            if( k < winEnd-2 )
                H(k+2,k-1) = zero;
        }
        else if( shiftStart > winBeg )
        {
            // The following is supposedly more reliable than
            // H(k,k-1) = -H(k,k-1) when v(1) and v(2) underflow
            // (cf. LAPACK's {s,d}lahqr)
            H(k,k-1) *= (one-tau0);
        }
        Real tau1 = tau0*v[1];
        if( numReflect == 3 )
        {
            Real tau2 = tau0*v[2];

            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd = H(k,j) + v[1]*H(k+1,j) + v[2]*H(k+2,j);
                H(k,  j) -= innerProd*tau0;
                H(k+1,j) -= innerProd*tau1;
                H(k+2,j) -= innerProd*tau2;
            }
           
            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+4,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd = H(j,k) + v[1]*H(j,k+1) + v[2]*H(j,k+2);
                H(j,k  ) -= innerProd*tau0;
                H(j,k+1) -= innerProd*tau1;
                H(j,k+2) -= innerProd*tau2;
            }

            if( ctrl.wantSchurVecs )
            {
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd = Z(j,k) + v[1]*Z(j,k+1) + v[2]*Z(j,k+2);
                    Z(j,k  ) -= innerProd*tau0;
                    Z(j,k+1) -= innerProd*tau1;
                    Z(j,k+2) -= innerProd*tau2;
                }
            }
        }
        else if( numReflect == 2 )
        {
            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd = H(k,j) + v[1]*H(k+1,j);
                H(k,  j) -= innerProd*tau0;
                H(k+1,j) -= innerProd*tau1;
            }

            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+3,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd = H(j,k) + v[1]*H(j,k+1);
                H(j,k  ) -= innerProd*tau0;
                H(j,k+1) -= innerProd*tau1;
            }

            if( ctrl.wantSchurVecs )
            {
                // Accumulate the Schur vectors
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd = Z(j,k) + v[1]*Z(j,k+1);
                    Z(j,k  ) -= innerProd*tau0;
                    Z(j,k+1) -= innerProd*tau1;
                }
            }
        }
    }
}

// Unfortunately, it seems to be the case that it is noticeably faster
// for this routine to manually inline the data access than to use the 
// (presumably inlined) Matrix::operator()(int,int) calls.
template<typename Real>
void SweepOpt
( Matrix<Real>& H,
  const Complex<Real>& shift0,
  const Complex<Real>& shift1,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    DEBUG_CSE
    const Real zero(0), one(1);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    Real* HBuf = H.Buffer();
    Real* ZBuf = Z.Buffer();
    const Int HLDim = H.LDim();
    const Int ZLDim = Z.LDim();

    const Int transformBeg = ( ctrl.fullTriangle ? 0 : winBeg ); 
    const Int transformEnd = ( ctrl.fullTriangle ? n : winEnd );

    vector<Real> v(3);
    auto subInd = IR(winBeg,winEnd);
    Int shiftStart = winBeg +
      ChooseStart( H(subInd,subInd), shift0, shift1, v );

    for( Int k=shiftStart; k<winEnd-1; ++k )
    {
        const Int numReflect = Min( 3, winEnd-k );
        if( k > shiftStart )
        {
            MemCopy( v.data(), &HBuf[k+(k-1)*HLDim], numReflect );
        }
        Real tau0 = lapack::Reflector( numReflect, v[0], &v[1], 1 );
        if( k > shiftStart )
        {
            HBuf[ k   +(k-1)*HLDim] = v[0];
            HBuf[(k+1)+(k-1)*HLDim] = zero;
            if( k < winEnd-2 )
                HBuf[(k+2)+(k-1)*HLDim] = zero;
        }
        else if( shiftStart > winBeg )
        {
            // The following is supposedly more reliable than
            // H(k,k-1) = -H(k,k-1) when v(1) and v(2) underflow
            // (cf. LAPACK's {s,d}lahqr)
            HBuf[k+(k-1)*HLDim] *= (one-tau0);
        }
        Real tau1 = tau0*v[1];
        if( numReflect == 3 )
        {
            Real tau2 = tau0*v[2];

            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd =
                       HBuf[ k   +j*HLDim] +
                  v[1]*HBuf[(k+1)+j*HLDim] +
                  v[2]*HBuf[(k+2)+j*HLDim];
                HBuf[ k   +j*HLDim] -= innerProd*tau0;
                HBuf[(k+1)+j*HLDim] -= innerProd*tau1;
                HBuf[(k+2)+j*HLDim] -= innerProd*tau2;
            }
           
            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+4,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd =
                       HBuf[j+ k   *HLDim] +
                  v[1]*HBuf[j+(k+1)*HLDim] +
                  v[2]*HBuf[j+(k+2)*HLDim];
                HBuf[j+ k   *HLDim] -= innerProd*tau0;
                HBuf[j+(k+1)*HLDim] -= innerProd*tau1;
                HBuf[j+(k+2)*HLDim] -= innerProd*tau2;
            }

            if( ctrl.wantSchurVecs )
            {
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd =
                           ZBuf[j+ k   *ZLDim] +
                      v[1]*ZBuf[j+(k+1)*ZLDim] +
                      v[2]*ZBuf[j+(k+2)*ZLDim];
                    ZBuf[j+ k   *ZLDim] -= innerProd*tau0;
                    ZBuf[j+(k+1)*ZLDim] -= innerProd*tau1;
                    ZBuf[j+(k+2)*ZLDim] -= innerProd*tau2;
                }
            }
        }
        else if( numReflect == 2 )
        {
            // Apply the Householder reflector from the left
            for( Int j=k; j<transformEnd; ++j )
            {
                Real innerProd = HBuf[k+j*HLDim] + v[1]*HBuf[(k+1)+j*HLDim];
                HBuf[ k   +j*HLDim] -= innerProd*tau0;
                HBuf[(k+1)+j*HLDim] -= innerProd*tau1;
            }

            // Apply the Householder reflector from the right
            const Int rightApplyEnd = Min(k+3,winEnd);
            for( Int j=transformBeg; j<rightApplyEnd; ++j )
            {
                Real innerProd = HBuf[j+k*HLDim] + v[1]*HBuf[j+(k+1)*HLDim];
                HBuf[j+ k   *HLDim] -= innerProd*tau0;
                HBuf[j+(k+1)*HLDim] -= innerProd*tau1;
            }

            if( ctrl.wantSchurVecs )
            {
                // Accumulate the Schur vectors
                for( Int j=0; j<nZ; ++j )
                {
                    Real innerProd = ZBuf[j+k*ZLDim] + v[1]*ZBuf[j+(k+1)*ZLDim];
                    ZBuf[j+ k   *ZLDim] -= innerProd*tau0;
                    ZBuf[j+(k+1)*ZLDim] -= innerProd*tau1;
                }
            }
        }
    }
}

} // namespace double_shift
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_DOUBLE_SHIFT_SWEEP_HPP
