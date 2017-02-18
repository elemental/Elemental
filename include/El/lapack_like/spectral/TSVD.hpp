/***
* 
* TSVD Based off of implementation by Andreas Noack. See: http://github.com/andreasnoack/TSVD.jl
*/
#ifndef EL_LAPACKLIKE_SPECTRAL_TSVD_HPP
#define EL_LAPACKLIKE_SPECTRAL_TSVD_HPP
#include <tuple>
#include <El/blas_like/level3.hpp>

namespace El {

template<typename Real>
Real CopySign(const Real& x, const Real& y) { return y >= Real(0) ? Abs(x) : -Abs(x); }

template<typename F>
struct BidiagonalMatrix {
 typedef Matrix <Base<F>> RealMatrix;
 
 BidiagonalMatrix(El::Int N) : mainDiag(N, 1), superDiag(N - 1, 1) {}
 
 RealMatrix mainDiag;
 RealMatrix superDiag;
}; //BidiagonalMatrix


struct BidiagInfo {
 Int numVecsReorth;
}; //struct BidiagInfo

template<typename F>
struct BidiagCtrl {
 bool reorthIn = false;
 Int minBlockSize=50;
 Base <F> reorthogTol = Pow(El::limits::Epsilon<El::Base<F>>(), El::Base<F>(.7));
 Base <F> convValTol = Pow(El::limits::Epsilon<El::Base<F>>(), El::Base<F>(.7));
 Base <F> convVecTol = Pow(El::limits::Epsilon<El::Base<F>>(), El::Base<F>(.5));
 bool partialProject = false;
}; //end struct BidiagCtrl

template<typename F>
struct ReorthogonalizationData {
 ReorthogonalizationData(El::Int N, F mu, F nu, F tau_): tau(tau_) {
  muList.reserve(N);
  nuList.reserve(N);
  muList.push_back(mu);
  muList.push_back(1);
  nuList.push_back(nu);
  nuList.push_back(1);
 }
 std::vector<El::Base<F>> muList, nuList;
 F tau;
};

//TODO: Make batch with Gemv
template<typename F>
Int Reorthogonalize(DistMatrix <F>& Q,
                    Int j,
                    std::vector<El::Base<F>>& termList,
                    const BidiagCtrl<F>& ctrl,
                    DistMatrix <F>& x,
                    Matrix <Base<F>>& diagList) {
 typedef Base<F> Real;
 const Real eps = limits::Epsilon<Real>();
 Int numVecsReorth = 0;
 termList.emplace_back(Real(1));
 for (Int i = 0; i < j; ++i) {
  auto qi = Q(ALL, IR(i));
  Axpy(-Dot(qi, x), qi, x);
  termList[i] = eps;
  numVecsReorth++;
 }
 Real alpha = Nrm2(x);
 diagList.Set(j, 0, alpha);
 Scale(Real(1) / alpha, x);
 return numVecsReorth;
}

template<typename F>
bool omegaRecurrenceV(const Int j,
                      const El::Base<F>& tau,
                      const BidiagonalMatrix <F>& B,
                      const std::vector<El::Base<F>>& muList,
                      const El::Base<F>& alpha,
                      const El::Base<F>& beta,
                      const BidiagCtrl<F>& ctrl,
                      std::vector<El::Base<F>>& nuList) {
 auto& superDiag = B.superDiag;
 auto& mainDiag = B.mainDiag;
 bool foundInaccurate = false; //run omega recurrence
 for (Int i = 0; i < j; ++i) {
  auto nu = superDiag.Get(i, 0) * muList[i + 1] +
            mainDiag.Get(i, 0) * muList[i] - beta * nuList[i];
  nu = ( nu + CopySign(tau, nu) ) / alpha;
  foundInaccurate |= ( Abs(nu) > ctrl.reorthogTol );
  nuList[i] = nu;
 }
 return foundInaccurate;
}


template<typename F>
bool omegaRecurrenceU(const Int j,
                      const El::Base<F>& tau,
                      const BidiagonalMatrix <F>& B,
                      const std::vector<El::Base<F>>& nuList,
                      const El::Base<F>& alpha,
                      const El::Base<F>& beta,
                      const BidiagCtrl<F>& ctrl,
                      std::vector<El::Base<F>>& muList) {
 auto& superDiag = B.superDiag;
 auto& mainDiag = B.mainDiag;
 bool foundInaccurate = false; //run omega recurrence
 for (Int i = 0; i <= j; ++i) {
  auto mu = mainDiag.Get(i, 0) * nuList[i] - alpha * muList[i];
  if ( i > 0 ) {
   mu += superDiag.Get(i - 1, 0) * nuList[i - 1];
  }
  mu = ( mu + CopySign(tau, mu) ) / beta;
  foundInaccurate |= ( Abs(mu) > ctrl.reorthogTol );
  muList[i] = mu;
 };
 return foundInaccurate;
}


template<typename F,
 typename ForwardOperator,
 typename AdjointOperator>
BidiagInfo
GolubKahan(
 El::Int iter,
 El::Int k,
 const ForwardOperator& A, const AdjointOperator& AAdj,
 BidiagonalMatrix <F>& B, El::DistMatrix<F>& U, El::DistMatrix<F>& V,
 ReorthogonalizationData<F>& reorthData, const El::BidiagCtrl<F>& ctrl) {
 
 typedef El::Base<F> Real;
 typedef El::DistMatrix<F> DM;
 const Real eps = El::limits::Epsilon<Real>();
 
 bool reorthB = ctrl.reorthIn;
 auto& nuList = reorthData.nuList;
 auto& muList = reorthData.muList;
 El::Int numVecsReorth = 0;
 //std::vector<El::Base<F>> maxMuList, maxNuList;
 //maxMuList.reserve(steps);
 //maxNuList.reserve(steps);
 
 DM v_j = V(ALL, IR(iter));
 DM u_j = U(ALL, IR(iter + 1));
 DM v_prev(v_j);
 DM u_prev(u_j);
 Real beta_j = 0;
 Real alpha = B.superDiag.Get(iter, 0);
 //auto absComp = [](const Real& a, const Real& b) { return Abs(a) < Abs(b); };
 for (Int j = iter + 1; j < k; ++j) {
 
  //The u step
  // apply the operator
  u_prev = u_j;
  u_j = A * v_prev;
  Axpy(-alpha, u_prev, u_j);
  beta_j = Nrm2(u_j);
  Scale(F(1.0)/beta_j, u_j);
  
  //reorthData.tau = Max(reorthData.tau, eps * ( alpha + beta ));
  ////compute omega recurrence
  //auto foundInaccurateU = omegaRecurrenceU(j, reorthData.tau, B, nuList, alpha, beta, ctrl, muList);
  //maxMuList.emplace_back(*std::max_element(muList.begin(), muList.end(), absComp));
  //muList.emplace_back(1);
  //if ( reorthB || foundInaccurateU ) {
  numVecsReorth += Reorthogonalize(U, j, muList, ctrl, u_j, B.superDiag);
  U(ALL, IR(j)) = u_j;
  //}
  
  v_prev = v_j;
  v_j = AAdj * u_j;
  Axpy(-beta_j, v_prev, v_j);
  Real alpha = Nrm2(v_j);
  Scale(F(1.0)/alpha, v_j);
  
  //reorthData.tau = Max(reorthData.tau, eps * ( alpha + beta ));
  //bool foundInaccurate;
  //Real maxElt;
  //auto foundInaccurateV = omegaRecurrenceV(j, reorthData.tau, B, muList, alpha, beta, ctrl, nuList);
  //maxNuList.emplace_back(*std::max_element(nuList.begin(), nuList.end(), absComp));
  ////Reorthogonalize if necessary.
  //if ( reorthB || foundInaccurateV ) {
  numVecsReorth += Reorthogonalize(V, j, nuList, ctrl, v_j, B.mainDiag);
  V(ALL, IR(j)) = v_j;
  //}
 }
 BidiagInfo info;
 info.numVecsReorth = numVecsReorth;
 return info;
}

template<typename F>
El::Matrix<F> BidiagonalSVD(const BidiagonalMatrix<F>& B,
                            Int n, Int numColsVT, Int numRowsU,
                            El::DistMatrix<F>& U, El::DistMatrix<F> VT) {
 auto s = B.mainDiag;
 auto superDiagCopy = B.superDiag;
 lapack::BidiagSVDQRAlg('U', n, numColsVT, numRowsU,
                        s.Buffer(), superDiagCopy.Buffer(),
                        VT.Buffer(), VT.LDim(),
                        U.Buffer(), U.LDim());
 return s;
}

template<typename F>
bool hasConverged(Int nVals,
                  Int newIndex,
                  const Matrix <Base<F>>& s,
                  const Matrix <Base<F>>& sOld,
                  const DistMatrix <F>& VTHat,
                  const DistMatrix <F>& UHat,
                  const El::Base<F>& beta,
                  const BidiagCtrl<F>& ctrl) {
 for (Int j = 0; j < nVals; ++j) {
  if ( Abs(s.Get(j, 0) - sOld.Get(j, 0)) > ctrl.convValTol ) {
   return false;
  }
  if ( beta * Abs(VTHat.Get(j, newIndex)) > ctrl.convVecTol ) {
   return false;
  }
  if ( beta * Abs(UHat.Get(newIndex, j)) > ctrl.convVecTol ) {
   return false;
  }
 }
 return true;
}

template< typename F,
          typename ForwardOperator,
          typename AdjointOperator>
ReorthogonalizationData<F> InitialIteration(const ForwardOperator& A, const AdjointOperator& AAdj,
                                            El::AbstractDistMatrix<F>& initialVec,
                                            DistMatrix<F>& UTilde, DistMatrix<F>& VTilde,
                                            BidiagonalMatrix<F>& B){
   typedef El::Base<F> Real;
   
   const Real eps = El::limits::Epsilon<Real>();
   
   Real initNorm = Nrm2(initialVec);
   Scale(Real(1)/initNorm, initialVec);
   UTilde(El::ALL, 0) = initialVec;
   DistMatrix<F> x = UTilde(El::ALL, IR(0));
   VTilde(El::ALL, 0) = AAdj * x;
   
   auto v = VTilde(El::ALL, 0);
  
   Real alpha = El::Nrm2(v); //record the norm, for reorth.
   El::Scale(Real(1)/alpha, v); //make v a unit vector
   
   UTilde(El::ALL, 1) = A * v;
   
   auto u0 = UTilde(El::ALL, IR(0));
   auto u1 = UTilde(El::ALL, IR(1));
   
   El::Axpy(-alpha, u0, u1);
   Real beta = El::Nrm2(u1);
   El::Scale(Real(1)/beta,u1);
   auto tau = eps*(alpha + beta);
   B.mainDiag.Set(0,0, alpha);
   B.superDiag.Set(0,0, beta);
   //1) Compute nu, a tolerance for reorthoganlization
   Real nu = 1 + tau/alpha; //think of nu = 1+\epsilon1
   Real mu = tau/beta;
   //first argument is maxIter which at the time of typing this was VTilde.Height and UTilde.Height().
   //If this changes revisit.
   ReorthogonalizationData<F> reorthData(VTilde.Height(), nu, mu, tau);
   
   return reorthData;
}
/**
* TSVD
*/
template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
inline El::Int
TSVD( 
      const ForwardOperator& A,
      const AdjointOperator& AAdj,
      El::Int nVals,
      El::AbstractDistMatrix<F>& U,
      El::AbstractDistMatrix<F>& S,
      El::AbstractDistMatrix<F>& V,
      El::AbstractDistMatrix<F>& initialVec, //should be of length m
      El::Int maxIter=El::Int(1000),
      const El::BidiagCtrl<F>& ctrl=El::BidiagCtrl<F>()) //requires InfiniteCharacteristic<F> && ForwardOperator::ScalarType is F
      {
        typedef El::DistMatrix<F> DM;
        typedef El::Matrix<El::Base<F>> RealMatrix;
        typedef El::Base<F> Real;
        
        const Real eps = El::limits::Epsilon<Real>();
        const El::Grid& g = initialVec.Grid();
        const El::Int m = A.Height();
        const El::Int n = A.Width();
       
        maxIter = El::Min(maxIter,El::Min(m,n));
 
        DM VTilde( n, maxIter, g);
        DM UTilde( m, maxIter, g);
 
        BidiagonalMatrix<F> B(maxIter+1);
      	 auto reorthData = InitialIteration(A, AAdj, initialVec, UTilde, VTilde, B);
        El::Display(UTilde, "UTilde");
        El::Display(VTilde, "VTilde");
        
        RealMatrix sOld( maxIter+1, 1);
        
        El::Int blockSize = 1;//El::Max(nVals, ctrl.minBlockSize);
        for(Int i = 0; i < maxIter; i+=blockSize){
            Int numSteps = Min( blockSize, maxIter-i);
            Int k = i+numSteps;
            //Write U_kA'V_k = B_k where k is i+numStep and B_k is diagonal.
            GolubKahan( i, k, A, AAdj, B, UTilde, VTilde, reorthData, ctrl);
            El::Display(UTilde(ALL, Range<Int>(0,k)), "UTilde");
            El::Display(VTilde(ALL, Range<Int>(0,k)), "VTilde");
         
            //C++17 this should be auto [UHat, S, VThat] = BidiagonalSVD(B,n,m);
            DM VTHat(k, nVals, g);
            DM UHat(nVals, k, g);
            auto s = BidiagonalSVD(B,k,nVals,nVals,UHat,VTHat);
            Display(UHat, "UHat");
            Display(VTHat, "VTHat");
            Real beta = B.superDiag.Get(k-1,0);
            if( hasConverged(nVals,  k, s, sOld,VTHat, UHat, beta, ctrl)){
                std::cout << "Apparently, we have converged" << std::endl;
                El::Gemm(El::NORMAL, El::NORMAL,  Real(1), UTilde,  UHat, Real(0), U);
                El::Gemm(El::NORMAL, El::ADJOINT, Real(1), VTilde, VTHat, Real(0), V);
                DistMatrix<F,STAR,STAR> S_STAR_STAR( S.Grid());
                S_STAR_STAR.LockedAttach( S.Height(), S.Width(), S.Grid(), S.ColAlign(), S.RowAlign(), s.Buffer(), s.LDim(), S.Root());
                Copy(S_STAR_STAR, S);
                return k;
            }
            sOld = s;
	           reorthData.tau = eps*s.Get(0,0);
        }
	return Int(-1);
}

/**
* TSVD don't provide initial guess.
*/
template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
inline Int 
TSVD( const ForwardOperator& A,
      const AdjointOperator& AAdj,
      Int nVals,
      AbstractDistMatrix<F>& U,
      AbstractDistMatrix<F>& S,
      AbstractDistMatrix<F>& V,
      Int maxIter=Int(1000),
      const BidiagCtrl<F>& ctrl=BidiagCtrl<F>()) {
    DistMatrix<F> initialVec(U.Grid());
    F center = rand();
    Uniform( initialVec, A.Height(), 1, center);
    return TSVD(A, AAdj, nVals, U,S,V, initialVec, maxIter, ctrl);
}

} //end namespace El
#endif //ifndef EL_LAPACKLIKE_SPECTRAL_TSVD
