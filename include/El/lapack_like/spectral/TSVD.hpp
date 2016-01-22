/***
* 
* TSVD Based off of implementation by Andreas Noack. See: http://github.com/andreasnoack/TSVD.jl
*/
#pragma once
#include <tuple>

namespace El {

template< typename Real>
Real CopySign(const Real& x, const Real& y){ return y >= Real(0) ? Abs(x): -Abs(x); }

struct BidiagInfo {
    Int numVecsReorth;
}; //struct BidiagInfo

template< typename F>
struct BidiagCtrl {
    bool reorthIn=false;
    Base<F> reorthogTol=Pow(limits::Epsilon<Base<F>>(), Base<F>(.7));
    Base<F> convValTol=Pow(limits::Epsilon<Base<F>>(), Base<F>(.7));
    Base<F> convVecTol=Pow(limits::Epsilon<Base<F>>(), Base<F>(.5));
    bool partialProject=false;
}; //end struct BidiagCtrl

//TODO: Make batch with Gemv
template< typename F>
Int reorth( DistMatrix<F>& Q,
            Int j,
            std::vector<Base<F>>& termList,
            const BidiagCtrl<F>& ctrl,
            DistMatrix<F>& x,
            Matrix<Base<F>>& diagList){
    typedef Base<F> Real;
    const Real eps = limits::Epsilon<Real>();
    Int numVecsReorth = 0;
    termList.emplace_back( Real( 1));
    for( Int i = 0; i < j; ++i){
        auto qi = Q(ALL, IR(i));
        Axpy( -Dot(qi, x), qi, x);
        termList[i] = eps;
        numVecsReorth++;
    }
    Real alpha = Nrm2( x);
    diagList.Set(j,0,alpha);
    Scale(Real(1)/alpha, x);
    auto qj = Q(ALL,IR(j));
    qj = x;
    return numVecsReorth;
}



template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
BidiagInfo
BidiagLanczos(
    Int m,
    Int n, 
    const ForwardOperator& A,
    const AdjointOperator& AAdj,
    Int steps,
    Base<F> tau,
    Matrix<Base<F>>& mainDiag,
    Matrix<Base<F>>& superDiag,
    DistMatrix<F>& U,
    DistMatrix<F>& V,
    //Base<F>& normA, //can be an estimate?
    std::vector<Base<F>>& muList,
    std::vector<Base<F>>& nuList,
    const BidiagCtrl<F>& ctrl){
    typedef Base<F> Real;
    typedef DistMatrix<F> DM;
    const Real eps = limits::Epsilon<Real>();

    bool reorthB = ctrl.reorthIn;
    Int numVecsReorth = 0;
    Int iter = mainDiag.size()-1;
    std::vector< Base<F>> maxMuList, maxNuList;
    maxMuList.reserve(steps);
    maxNuList.reserve(steps);

    DM v = V( ALL, IR(iter));
    DM u = U( ALL, IR(iter+1));
    DM vOld( v);
    DM uOld( u);
    Real beta = superDiag[iter];
    for( Int j = iter+1; j <= iter+steps; ++j)
    {
        vOld = v;
        v = u;
        AAdj( v);
        Axpy( -beta, vOld, v);
        Real alpha = Nrm2(v);
        //run omega recurrence
        bool foundInaccurate = false;
        for(Int i = 0; i < j; ++i)
        {
            Real nu = superDiag.Get(i,0)*muList[i+1] + 
                      mainDiag.Get(i,0)*muList[i] - beta*nuList[i];
            nu = (nu + CopySign(tau, nu))/alpha;
            foundInaccurate |= (Abs(nu) > ctrl.reorthogTol); 
            nuList[i] = nu;
        }
        Real maxElt = *std::max_element( nuList.begin(), nuList.end(), Abs);
        maxNuList.emplace_back( maxElt);
        if( reorthB || foundInaccurate ){
            numVecsReorth += reorth( V, j, nuList, ctrl, v, mainDiag);
        }

        //The u step
        uOld = u;
        u = v;
        // apply the operator
        A(u);
        Axpy(-alpha, uOld, u);
        beta = Nrm2( u);

        //compute omega recurrence
        foundInaccurate=false;
        for( Int i = 0; i <= j; ++i){
            Real mu = mainDiag.Get(i,0)*nuList[i] - alpha*muList[i];
            if (i > 0){
                mu += superDiag.Get(i-1, 0)*nuList[i-1];
            }
            mu = (mu + CopySign(tau, mu))/beta;
            foundInaccurate |= (Abs(mu) > ctrl.reorthogTol); 
            muList[ i] = mu;
        }
        maxElt = *std::max_element( muList.begin(), muList.end(), Abs);
        maxMuList.emplace_back( maxElt);
        muList.emplace_back( 1);
        if( reorthB || foundInaccurate){
            numVecsReorth += reorth( U, j, muList, ctrl, u, superDiag); 
        }
    }
    BidiagInfo info;
    info.numVecsReorth = numVecsReorth;
    return info;
}


/**
* TSVD
*/
template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
inline Int 
tsvd( 
      Int m,
      Int n, 
      const ForwardOperator& A,
      const AdjointOperator& AAdj,
      Int nVals, 
      DistMatrix<F>& initialVec, 
      Int maxIter=Int(1000),
      const BidiagCtrl<F>& ctrl=BidiagCtrl<F>()) //requires InfiniteCharacteristic<F> && ForwardOperator::ScalarType is F
      {
        typedef DistMatrix<F> DM;
        typedef Matrix<Base<F>> RealMatrix;
        typedef Base<F> Real;
        maxIter = Min(maxIter,Min(m,n));
        const Grid& g = initialVec.Grid();
        //1) Compute nu, a tolerance for reorthoganlization
        auto tau = InfinityNorm(A); //TODO: Impl this?
        
        Real initNorm = Nrm2(initialVec);
        Scale(Real(1)/initNorm, initialVec);       
        DM V( n, maxIter, g);
        auto v = V(ALL, 0);
        v = initialVec;
        AAdj( v); //v = A'initialVec;
        Real alpha = Nrm2(v); //record the norm, for reorth.
        Scale(Real(1)/alpha, v); //make v a unit vector
        Real nu = 1 + tau/alpha; //think of nu = 1+\epsilon
        
        DM U(m, maxIter, g);
        auto u0 = U(ALL, 0);
        auto u1 = U(ALL, 1);
        u0 = initialVec;
        u1 = v;
        A( u1);
        Axpy(-alpha, u0, u1);
        Real beta = Nrm2(u1);
        Scale(Real(1)/beta,u1);
        Real mu = tau/beta;


        RealMatrix mainDiag( maxIter+1, 1); //TODO: make this a Matrix<Real>
        RealMatrix superDiag( maxIter, 1, g); //TODO: make this a Matrix<Real>
        //TODO: Implement the rest.
        std::vector<Base<F>> muList, nuList;
        muList.reserve( maxIter);
        nuList.reserve( maxIter);
        muList.push_back( mu);
        muList.push_back( 1);
        nuList.push_back( 1);
       
        RealMatrix sOld( maxIter+1, 1);
        Int blockSize = Max(nVals, 50);
        for(Int i = 0; i < maxIter; i+=blockSize){
            Int numSteps = Min(blockSize,maxIter-i); 
            BidiagLanczos
            ( m, n, A, AAdj, numSteps, 
              tau, mainDiag, superDiag, U, V, 
              muList, nuList, ctrl);
            auto s = mainDiag;
            auto superDiagCopy = superDiag;
            RealMatrix VT(i+numSteps,nVals);
            RealMatrix U(nVals,i+numSteps);
            lapack::BidiagQRAlg
            (
             'U', i+numSteps, nVals, nVals,
             s.Buffer(), superDiagCopy.Buffer(), 
             VT.Buffer(), VT.LDim(),
             U.Buffer(), U.LDim());

            Real beta = superDiag.Get(i+numSteps-1,0);
            bool converged=true;
            for(Int j = 0; j < nVals; ++j){
                if( Abs(s.Get(j,0) - sOld.Get(j,0)) > ctrl.convValTol){
                    converged = false;
                }
                if( beta*Abs(VT.Get(j, i+numSteps)) > ctrl.convVecTol){
                    converged = false;
                }
                if( beta*Abs(U.Get(i+numSteps, j)) > ctrl.convVecTol){
                    converged=false;
                }
            }
            sOld = s;
        }
}
} //end namespace El
