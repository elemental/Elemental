#pragma once
#include <tuple>

namespace El {

template< typename F, 
          typename ForwardOperator,
          typename AdjointOperator>
std::tuple< bool, Base<F>, Base<F> >
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
    Base<F> tolError){
    typedef Base<F> Real;
    typename DistMatrix<F> DM;

    bool reorthB = reorthIn;
    Int nReorth = 0;
    Int nReorthVec = 0;
    muList.reserve( steps);
    nuList.reserve( steps);
    Int iter = alphaList.size()-1;
    DM u = U(ALL, IR(iter+1));
    DM v = V(ALL, IR(iter));
    Real beta = betaList[iter];

    
    DM vOld( v);
    for( Int j = iter+1; j <= iter+steps; ++j)
    {
        vOld = v;
        v = u;
        AAdj( v);
        Axpy( -beta, vOld, v);
        Real alpha = Nrm2(v);
        //run omega recurrence
        for(Int i = 0; i < j; ++i)
        {
            
        }

    }

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
        Int cc = Max((Int)5, nVals);
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

}
} //end namespace El



