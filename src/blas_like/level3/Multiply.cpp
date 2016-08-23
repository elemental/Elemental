/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level3.hpp>

namespace El {

namespace {
/**
 * MultiplyCSR specialization where the CSR matrix happens to have all nonzeros = 1.
 */
template<typename T>
void MultiplyCSR
( Orientation orientation,
  Int m, Int n,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   x,
  T beta,
        T*   y )
{
    DEBUG_CSE
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            T sum = 0;
            const Int eStart = rowOffsets[i];
            const Int eStop = rowOffsets[i+1];
            for( Int e=eStart; e<eStop; ++e )
                sum += x[colIndices[e]];         
            y[i] = alpha*sum + beta*y[i];
        }
    }
    else
    {
        for( Int j=0; j<n; ++j ){
            y[j] *= beta;
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    y[colIndices[e]] += alpha*x[i];         
            }
	}
    }
}

template<typename T,typename=DisableIf<IsBlasScalar<T>>>
void MultiplyCSR
( Orientation orientation,
  Int m, Int n,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   x,
  T beta,
        T*   y )
{
    DEBUG_CSE
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            T sum = 0;
            const Int eStart = rowOffsets[i];
            const Int eStop = rowOffsets[i+1];
            for( Int e=eStart; e<eStop; ++e )
                sum += values[e]*x[colIndices[e]];         
            y[i] = alpha*sum + beta*y[i];
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int j=0; j<n; ++j )
            y[j] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
            }
        }
    }
}

template<typename T,typename=EnableIf<IsBlasScalar<T>>,typename=void>
void MultiplyCSR
( Orientation orientation,
  Int m, Int n,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   x,
  T beta,
        T*   y )
{
    DEBUG_CSE
#if defined(EL_HAVE_MKL) && !defined(EL_DISABLE_MKL_CSRMV)
    char matDescrA[6];
    matDescrA[0] = 'G';
    matDescrA[3] = 'C';
    mkl::csrmv
    ( orientation, m, n, alpha, matDescrA, 
      values, colIndices, rowOffsets, rowOffsets+1, x, beta, y );
#else
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            T sum = 0;
            const Int eStart = rowOffsets[i];
            const Int eStop = rowOffsets[i+1];
            for( Int e=eStart; e<eStop; ++e )
                sum += values[e]*x[colIndices[e]];         
            y[i] = alpha*sum + beta*y[i];
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int j=0; j<n; ++j )
            y[j] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
            }
        }
    }
#endif
}

template<typename T>
void MultiplyCSR
( Orientation orientation,
  Int m, Int n, Int numRHS,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   X, Int ldX,
  T beta,
        T*   Y, Int ldY )
{
    DEBUG_CSE
    if( numRHS == 1 )
    {
        MultiplyCSR
        ( orientation, m, n, alpha, 
          rowOffsets, colIndices, values, X, beta, Y );
        return;
    }

    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            for( Int k=0; k<numRHS; ++k )
            {
                T sum = 0;
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    sum += values[e]*X[colIndices[e]+k*ldX];
                Y[i+k*ldY] = alpha*sum + beta*Y[i+k*ldY];
            }
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int k=0; k<numRHS; ++k )
            for( Int j=0; j<n; ++j )
                Y[j+k*ldY] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i+k*ldX];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i+k*ldX];         
                }
            }
        }
    }
}

template<typename T>
void MultiplyCSR
( Orientation orientation,
  Int m, Int n, Int numRHS,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   X, Int ldX,
  T beta,
        T*   Y, Int ldY )
{
    DEBUG_CSE
    if( numRHS == 1 )
    {
        MultiplyCSR
        ( orientation, m, n, alpha, 
          rowOffsets, colIndices, X, beta, Y );
        return;
    }

    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            for( Int k=0; k<numRHS; ++k )
            {
                T sum = 0;
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    sum += X[colIndices[e]+k*ldX];
                Y[i+k*ldY] = alpha*sum + beta*Y[i+k*ldY];
            }
        }
    }
    else
    {
        for( Int k=0; k<numRHS; ++k )
	{
            for( Int j=0; j<n; ++j )
	    {
                Y[j+k*ldY] *= beta;
            	for( Int i=0; i<m; ++i )
            	{
            	    const Int eStart = rowOffsets[i];
            	    const Int eStop = rowOffsets[i+1];
            	    for( Int e=eStart; e<eStop; ++e )
            	    {
            	        for( Int k=0; k<numRHS; ++k ){
            	            Y[colIndices[e]+k*ldY] += alpha*X[i+k*ldX]; 
			}	    
            	    }
            	}
            }
        }
    }
}

template<typename T>
void MultiplyCSRInterX
( Orientation orientation,
  Int m, Int n, Int numRHS,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   X,
  T beta,
        T*   Y, Int ldY )
{
    DEBUG_CSE
    if( numRHS == 1 )
    {
        MultiplyCSR
        ( orientation, m, n, alpha, 
          rowOffsets, colIndices, values, X, beta, Y );
        return;
    }

    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            for( Int k=0; k<numRHS; ++k )
            {
                T sum = 0;
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    sum += values[e]*X[colIndices[e]*numRHS+k];
                Y[i+k*ldY] = alpha*sum + beta*Y[i+k*ldY];
            }
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int k=0; k<numRHS; ++k )
            for( Int j=0; j<n; ++j )
                Y[j+k*ldY] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i*numRHS+k];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i*numRHS+k];         
                }
            }
        }
    }
}

template<typename T>
void MultiplyCSRInterY
( Orientation orientation,
  Int m, Int n, Int numRHS,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   X, Int ldX,
  T beta,
        T*   Y )
{
    DEBUG_CSE
    if( numRHS == 1 )
    {
        MultiplyCSR
        ( orientation, m, n, alpha, 
          rowOffsets, colIndices, values, X, beta, Y );
        return;
    }

    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            for( Int k=0; k<numRHS; ++k )
            {
                T sum = 0;
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    sum += values[e]*X[colIndices[e]+k*ldX];
                Y[i*numRHS+k] = alpha*sum + beta*Y[i*numRHS+k];
            }
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int k=0; k<numRHS; ++k )
            for( Int j=0; j<n; ++j )
                Y[j*numRHS+k] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i+k*ldX];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i+k*ldX];         
                }
            }
        }
    }
}

template<typename T>
void MultiplyCSRInter
( Orientation orientation,
  Int m, Int n, Int numRHS,
  T alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const T*   values,
  const T*   X,
  T beta,
        T*   Y )
{
    DEBUG_CSE
    if( numRHS == 1 )
    {
        MultiplyCSR
        ( orientation, m, n, alpha, 
          rowOffsets, colIndices, values, X, beta, Y );
        return;
    }

    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            for( Int k=0; k<numRHS; ++k )
            {
                T sum = 0;
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                    sum += values[e]*X[colIndices[e]*numRHS+k];
                Y[i*numRHS+k] = alpha*sum + beta*Y[i*numRHS+k];
            }
        }
    }
    else
    {
        const bool conj = ( orientation == ADJOINT );
        for( Int j=0; j<n; ++j )
            for( Int k=0; k<numRHS; ++k )
                Y[j*numRHS+k] *= beta;
        if( conj )
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i*numRHS+k];
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                const Int eStart = rowOffsets[i];
                const Int eStop = rowOffsets[i+1];
                for( Int e=eStart; e<eStop; ++e )
                {
                    T prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i*numRHS+k]; 
                }
            }
        }
    }
}

} // anonymous namespace

template<typename T>
void Multiply
( Orientation orientation,
  T alpha, const SparseMatrix<T>& A, const Matrix<T>& X,
  T beta,                                  Matrix<T>& Y )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( X.Width() != Y.Width() )
          LogicError("X and Y must have the same width");
    )
    MultiplyCSR
    ( orientation, A.Height(), A.Width(), X.Width(),
      alpha, A.LockedOffsetBuffer(),
             A.LockedTargetBuffer(),
             A.LockedValueBuffer(),
             X.LockedBuffer(), X.LDim(),
      beta,  Y.Buffer(),       Y.LDim() );
}

template<typename T>
void Multiply
( Orientation orientation, 
  T alpha, const Graph& A, const Matrix<T>& X,
  T beta,                        Matrix<T>& Y )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( X.Width() != Y.Width() )
          LogicError("X and Y must have the same width");
    )
    MultiplyCSR
    ( orientation, A.NumSources(), A.NumTargets(), X.Width(), 
      alpha, A.LockedOffsetBuffer(), 
             A.LockedTargetBuffer(), 
             X.LockedBuffer(), X.LDim(),
      beta,  Y.Buffer(), Y.LDim());
}


template<typename T>
void Multiply
( Orientation orientation, 
        T alpha, 
  const DistSparseMatrix<T>& A,
  const DistMultiVec<T>& X,
        T beta,
        DistMultiVec<T>& Y )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( X.Width() != Y.Width() )
          LogicError("X and Y must have the same width");
      if( !mpi::Congruent( A.Comm(), X.Comm() ) || 
          !mpi::Congruent( X.Comm(), Y.Comm() ) )
          LogicError("Communicators did not match");
    )

    const bool time = false;

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    // TODO: Use sequential implementation if commSize = 1?

    Timer totalTimer, timer;
    if( time && commRank == 0 )
        totalTimer.Start();

    // Y := beta Y
    Y *= beta;

    A.InitializeMultMeta();
    const auto& meta = A.LockedDistGraph().multMeta;
    // Convert the sizes and offsets to be compatible with the current width
    const Int b = X.Width();
    vector<int> recvSizes=meta.recvSizes,
                recvOffs=meta.recvOffs,
                sendSizes=meta.sendSizes,
                sendOffs=meta.sendOffs;
    for( int q=0; q<commSize; ++q )
    {
        recvSizes[q] *= b;    
        recvOffs[q] *= b;
        sendSizes[q] *= b;
        sendOffs[q] *= b;
    }

    if( orientation == NORMAL )
    {
        if( A.Height() != Y.Height() )
            LogicError("A and Y must have the same height");
        if( A.Width() != X.Height() )
            LogicError("The width of A must match the height of X");

        // Pack the send values
        const Int numSendInds = meta.sendInds.size();
        const Int firstLocalRow = X.FirstLocalRow();
        vector<T> sendVals;
        FastResize( sendVals, numSendInds*b );
        const T* XBuffer = X.LockedMatrix().LockedBuffer();
        const Int ldX = X.LockedMatrix().LDim();
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                sendVals[s*b+t] = XBuffer[iLoc+t*ldX];
        }

        // Now send them
        vector<T> recvVals( meta.numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), sendSizes.data(), sendOffs.data(),
          recvVals.data(), recvSizes.data(), recvOffs.data(), comm );
     
        // Perform the local multiply-accumulate, y := alpha A x + y
        if( time && commRank == 0 )
            timer.Start();
        MultiplyCSRInterX
        ( NORMAL, A.LocalHeight(), meta.numRecvInds, b,
          alpha, A.LockedOffsetBuffer(), 
                 meta.colOffs.data(),
                 A.LockedValueBuffer(),
                 recvVals.data(), 
          T(1),  Y.Matrix().Buffer(), Y.Matrix().LDim() );
        if( time && commRank == 0 )
            Output("  MultiplyCSRInterX time: ",timer.Stop());
    }
    else
    {
        if( A.Width() != Y.Height() )
            LogicError("The width of A must match the height of Y");
        if( A.Height() != X.Height() )
            LogicError("The height of A must match the height of X");

        // Form and pack the updates to Y
        if( time && commRank == 0 )
            timer.Start();
        vector<T> sendVals( meta.numRecvInds*b, 0 );
        MultiplyCSRInterY
        ( orientation, A.LocalHeight(), meta.numRecvInds, b,
          alpha, A.LockedOffsetBuffer(),
                 meta.colOffs.data(),
                 A.LockedValueBuffer(),
                 X.LockedMatrix().LockedBuffer(), X.LockedMatrix().LDim(),
          T(1),  sendVals.data() );
        if( time && commRank == 0 )
            Output("  MultiplyCSRInterY time: ",timer.Stop());

        // Inject the updates to Y into the network
        const Int numRecvInds = meta.sendInds.size();
        vector<T> recvVals;
        FastResize( recvVals, numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), recvSizes.data(), recvOffs.data(),
          recvVals.data(), sendSizes.data(), sendOffs.data(), comm );
     
        // Accumulate the received indices onto Y
        const Int firstLocalRow = Y.FirstLocalRow();
        T* YBuffer = Y.Matrix().Buffer(); 
        const Int ldY = Y.Matrix().LDim();
        for( Int s=0; s<numRecvInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                YBuffer[iLoc+t*ldY] += recvVals[s*b+t];
        }
    }
    if( time && commRank == 0 )
        Output("Multiply total time: ",totalTimer.Stop());
}

#define PROTO(T) \
    template void Multiply \
    ( Orientation orientation, \
            T alpha, \
      const SparseMatrix<T>& A, \
      const Matrix<T>& X, \
            T beta, \
            Matrix<T>& Y ); \
    template void Multiply \
    ( Orientation orientation, \
            T alpha, \
      const Graph& A, \
      const Matrix<T>& X, \
            T beta, \
            Matrix<T>& Y ); \
    template void Multiply \
    ( Orientation orientation, \
            T alpha, \
      const DistSparseMatrix<T>& A, \
      const DistMultiVec<T>& X, \
            T beta, \
            DistMultiVec<T>& Y );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
