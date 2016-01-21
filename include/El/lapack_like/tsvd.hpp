#pragma once
#include <tuple>

namespace El {

template< typename Real>
Real CopySign(const Real& x, const Real& y){ return y >= Real(0) ? Abs(x): -Abs(x); }

template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
void
BidiagLanczos(
    Int m,
    Int n, 
    const ForwardOperator& A,
    const AdjointOperator& AAdj,
    Int steps,
    Base<F> tau,
    std::vector<Base<F>>& alphaList,
    std::vector<Base<F>>& betaList,
    DistMatrix<F>& U,
    DistMatrix<F>& V,
    //Base<F>& normA, //can be an estimate?
    std::vector<Base<F>>& muList,
    std::vector<Base<F>>& nuList,
    bool reorthIn,
    Base<F> tolError,
    bool partialProject=false){
    typedef Base<F> Real;
    typename DistMatrix<F> DM;
    const Real eps = limits::Epsilon<Real>();

    bool reorthB = reorthIn;
    Int nReorth = 0;
    Int nReorthVec = 0;
    muList.reserve( steps);
    nuList.reserve( steps);
    Int iter = alphaList.size()-1;
    DM u = U(ALL, IR(iter+1));
    DM v = V(ALL, IR(iter));
    Real beta = betaList[iter];

    
    std::vector< Int> reorth_nu;
    if( partialProject){ reorth_nu.reserve( steps); }
    DM vOld( v);
    for( Int j = iter+1; j <= iter+steps; ++j)
    {
        vOld = v;
        v = u;
        AAdj( v);
        Axpy( -beta, vOld, v);
        Real alpha = Nrm2(v);
        //run omega recurrence
        reorth_nu.clear();
        for(Int i = 0; i < j; ++i)
        {
            Real nu = betaList[i]*muList[i+1] + 
                      alphaList[i]*muList[i] - beta*nuList[i];
            nu = (nu + CopySign(tau, nu))/alpha;
            if( partialProject && Abs(v) > tolError) { 
                reorth_nu.emplace_back( i);
            }
            nuList[i] = nu;
        }
        Real maxElt = *std::max_element( nuList.begin(), nuList.end(), Abs);
        maxNuList.emplace_back( maxElt);
    }
    nuList.emplace_back( Real(1));
    if (reorth_b || reorthNu.size()>0){
        if( partialProject){
            for( Int i: reorth_nu){
                auto vi = V(ALL, IR(i));
                Axpy( -Dot(vi, v), vi, v);
                nuList[i] = 2*eps;
                nReorthVecs++;
            }
        } else { //full reorthogalization
             for( Int i = 0; i < j; ++i){
                auto vi = V(ALL, IR(i));
                Axpy( -Dot(vi, v), vi, v);
                nuList[i] = eps;
                nReorthVecs++;
            }
        }
        alpha = Nrm2( v);
        reorth_b = !reorth_b;
    }
    alphaList.emplace_back(alpha);
    Scale(Real(1)/alpha, v);
    auto vj = V(ALL,IR(j));
    vj = v;

    //The u step
    uOld = u;
    u = v;
    // apply the operator
    A(u);
    Axpy(-alpha, uOld, u);
    beta = Nrm2( u);

    //compute omega recurrence
    std::vector< Int> reorth_mu;
    if( partialProject){ reorth_mu.reserve( steps); }
    for( Int i = 0; i <= j; ++i){
        Real mu = alphaList[i]*nuList[i] - alpha*muList[i];
        if (i > 0){
            mu += betaList[i-1]*nuList[i-1];
        }
        mu = (mu + CopySign(tau, mu))/beta;
        if( partialProject && Abs(mu) > tolError){
            reorth_mu.emplace_back( i);
        }
        muList[ i] = mu;
    }
    Real maxElt = *std::max_element( muList.begin(), muList.end(), Abs);
    maxMuList.emplace_back( maxElt);
    muList.emplace_back( 1);
    
    if (reorth_b || reorthMu.size()>0){
        if( partialProject){
            for( Int i: reorth_mu){
                auto ui = U(ALL, IR(i));
                Axpy( -Dot(ui, u), ui, u);
                muList[i] = 2*eps;
                nReorthVecs++;
            }
        } else { //full reorthogalization
             for( Int i = 0; i < j; ++i){
                auto ui = U(ALL, IR(i));
                Axpy( -Dot(ui, u), ui, u);
                muList[i] = eps;
                nReorthVecs++;
            }
        }
        alpha = Nrm2( u);
        reorth_b = !reorth_b;
        nReorth++;
    }
    betaList.emplace_back(beta);
    Scale(Real(1)/beta, u);
    auto uj = U(ALL,IR(j));
    uj = u;
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
      Int maxIter, 
      DistMatrix<F>& initialVec, 
      Base<F> tolConv=Pow(limits::Epsilon<Base<F>>(), Base<F>(.7)),
      Base<F> tolError=Pow(limits::Epsilon<Base<F>>(), Base<F>(.7)),
      ) //requires InfiniteCharacteristic<F> && ForwardOperator::ScalarType is F
      {
        typedef DistMatrix<F> DM;
        typedef Base<F> Real;
        const Grid& g = initialVec.Grid();
        auto tau = InfinityNorm(A); //TODO: Impl this?
        DM z( g);
        Zeros(z, 1,1);
        Int steps = 5;
        Real initNorm = Nrm2(initial_vec);
        DM v( initialVec);
        AAdj( v);
        Scale(Real(1)/initNorm, v);       
        Real alpha = Nrm2(v);
        Scale(Real(1)/alpha, v);
        DM V( n, 1, g);
        Fill( V, Real(1));
        //TODO: Make this the right size upfront?
        DM alphaList( 1, 1, g); //TODO: make this a Matrix<Real>
        
        Real nu = 1 + tau/alpha; //think of nu = 1+\epsilon
        DM u( v);
        A( u);
        DM uOld( initVec); //TODO: similar?
        Scale(Real(1)/initNorm, initVec);       
        Axpy(-alpha, uOld, u);
        Real beta = FrobeniusNorm(u);
        Scale(Real(1)/beta,u);

        DM U(g);
        HCat(uOld, u, U);
        //TODO: Make this the right size upfront
        DM betaList( 1, 1, g); //TODO: make this a Matrix<Real>
        Real mu = tau/beta;
        std::vector< Real> maxMuList; 
        std::vector< Real> maxNuList; 
        maxMuList.reserve(maxIter);
        maxNuList.reserve(maxIter);
        ///TODO: Implement the rest.
        BidiagLanczos( m, n, A, AAdj, numSteps, 
                       tau, alphaList, betaList, U, V, 
                       muList, nuList, 
                       reorthIn, tolError);
        for( 


}
} //end namespace El



