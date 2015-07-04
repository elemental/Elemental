/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace {

template<typename Ring>
void MultiplyCSR
( Orientation orientation, Int m, Int n,
  Ring alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Ring*   values,
  const Ring*   x,
  Ring beta,
        Ring*   y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSR"))
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
            Ring sum = 0;
            for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
        }
        else
        {
            for( Int i=0; i<m; ++i )
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
        }
    }
#endif
}

template<>
void MultiplyCSR<Int>
( Orientation orientation, Int m, Int n,
  Int alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Int*   values,
  const Int*   x,
  Int beta,
        Int*   y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSR"))
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            Int sum = 0;
            for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
        }
        else
        {
            for( Int i=0; i<m; ++i )
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
        }
    }
}

#ifdef EL_HAVE_QUAD
template<>
void MultiplyCSR<Quad>
( Orientation orientation, Int m, Int n,
  Quad alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Quad*   values,
  const Quad*   x,
  Quad beta,
        Quad*   y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSR"))
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            Quad sum = 0;
            for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
        }
        else
        {
            for( Int i=0; i<m; ++i )
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
        }
    }
}

template<>
void MultiplyCSR<Complex<Quad>>
( Orientation orientation, Int m, Int n,
  Complex<Quad> alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Complex<Quad>*   values,
  const Complex<Quad>*   x,
  Complex<Quad> beta,
        Complex<Quad>*   y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSR"))
    if( orientation == NORMAL )
    {
        for( Int i=0; i<m; ++i )
        {
            Complex<Quad> sum = 0;
            for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*Conj(values[e])*x[i];         
        }
        else
        {
            for( Int i=0; i<m; ++i )
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                    y[colIndices[e]] += alpha*values[e]*x[i];         
        }
    }
}
#endif // ifdef EL_HAVE_QUAD

template<typename Ring>
void MultiplyCSR
( Orientation orientation, Int m, Int n, Int numRHS,
  Ring alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Ring*   values,
  const Ring*   X, Int ldX,
  Ring beta,
        Ring*   Y, Int ldY )
{
    DEBUG_ONLY(CSE cse("MultiplyCSR"))
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
                Ring sum = 0;
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i+k*ldX];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i+k*ldX];         
                }
            }
        }
    }
}

template<typename Ring>
void MultiplyCSRInterX
( Orientation orientation, Int m, Int n, Int numRHS,
  Ring alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Ring*   values,
  const Ring*   X,
  Ring beta,
        Ring*   Y, Int ldY )
{
    DEBUG_ONLY(CSE cse("MultiplyCSRInterX"))
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
                Ring sum = 0;
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i*numRHS+k];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]+k*ldY] += prod*X[i*numRHS+k];         
                }
            }
        }
    }
}

template<typename Ring>
void MultiplyCSRInterY
( Orientation orientation, Int m, Int n, Int numRHS,
  Ring alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Ring*   values,
  const Ring*   X, Int ldX,
  Ring beta,
        Ring*   Y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSRInterY"))
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
                Ring sum = 0;
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i+k*ldX];         
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i+k*ldX];         
                }
            }
        }
    }
}

template<typename Ring>
void MultiplyCSRInter
( Orientation orientation, Int m, Int n, Int numRHS,
  Ring alpha,
  const Int* rowOffsets,
  const Int* colIndices,
  const Ring*   values,
  const Ring*   X,
  Ring beta,
        Ring*   Y )
{
    DEBUG_ONLY(CSE cse("MultiplyCSRInter"))
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
                Ring sum = 0;
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
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
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*Conj(values[e]);
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i*numRHS+k];
                }
            }
        }
        else
        {
            for( Int i=0; i<m; ++i )
            {
                for( Int e=rowOffsets[i]; e<rowOffsets[i+1]; ++e )
                {
                    Ring prod = alpha*values[e];
                    for( Int k=0; k<numRHS; ++k )
                        Y[colIndices[e]*numRHS+k] += prod*X[i*numRHS+k]; 
                }
            }
        }
    }
}

} // anonymous namespace

template<typename Ring>
void Multiply
( Orientation orientation, 
  Ring alpha, const SparseMatrix<Ring>& A, 
              const Matrix<Ring>& X,
  Ring beta,        Matrix<Ring>& Y )
{
    DEBUG_ONLY(
      CSE cse("Multiply");
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

template<typename Ring>
void Multiply
( Orientation orientation, 
  Ring alpha, const DistSparseMatrix<Ring>& A, 
              const DistMultiVec<Ring>& X,
  Ring beta,        DistMultiVec<Ring>& Y )
{
    DEBUG_ONLY(
      CSE cse("Multiply");
      if( X.Width() != Y.Width() )
          LogicError("X and Y must have the same width");
      if( !mpi::Congruent( A.Comm(), X.Comm() ) || 
          !mpi::Congruent( X.Comm(), Y.Comm() ) )
          LogicError("Communicators did not match");
    )
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    // TODO: Use sequential implementation if commSize = 1?

    // Y := beta Y
    Y *= beta;

    A.InitializeMultMeta();
    const auto& meta = A.multMeta;

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
        vector<Ring> sendVals( numSendInds*b );
        const Ring* XBuffer = X.LockedMatrix().LockedBuffer();
        const Int ldX = X.LockedMatrix().LDim();
        for( Int s=0; s<numSendInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                sendVals[s*b+t] = XBuffer[iLoc+t*ldX];
        }

        // Now send them
        vector<Ring> recvVals( meta.numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), sendSizes.data(), sendOffs.data(),
          recvVals.data(), recvSizes.data(), recvOffs.data(), comm );
     
        // Perform the local multiply-accumulate, y := alpha A x + y
        MultiplyCSRInterX
        ( NORMAL, A.LocalHeight(), A.Width(), X.Width(),
          alpha, A.LockedOffsetBuffer(), 
                 meta.colOffs.data(),
                 A.LockedValueBuffer(),
                 recvVals.data(), 
          Ring(1), Y.Matrix().Buffer(), Y.Matrix().LDim() );
    }
    else
    {
        if( A.Width() != Y.Height() )
            LogicError("The width of A must match the height of Y");
        if( A.Height() != X.Height() )
            LogicError("The height of A must match the height of X");

        // Form and pack the updates to Y
        vector<Ring> sendVals( meta.numRecvInds*b, 0 );
        MultiplyCSRInterY
        ( orientation, A.LocalHeight(), A.Width(), X.Width(),
          alpha, A.LockedOffsetBuffer(),
                 meta.colOffs.data(),
                 A.LockedValueBuffer(),
                 X.LockedMatrix().LockedBuffer(), X.LockedMatrix().LDim(),
          Ring(1), sendVals.data() );

        // Inject the updates to Y into the network
        const Int numRecvInds = meta.sendInds.size();
        vector<Ring> recvVals( numRecvInds*b );
        mpi::AllToAll
        ( sendVals.data(), recvSizes.data(), recvOffs.data(),
          recvVals.data(), sendSizes.data(), sendOffs.data(), comm );
     
        // Accumulate the received indices onto Y
        const Int firstLocalRow = Y.FirstLocalRow();
        Ring* YBuffer = Y.Matrix().Buffer(); 
        const Int ldY = Y.Matrix().LDim();
        for( Int s=0; s<numRecvInds; ++s )
        {
            const Int i = meta.sendInds[s];
            const Int iLoc = i - firstLocalRow;
            for( Int t=0; t<b; ++t )
                YBuffer[iLoc+t*ldY] += recvVals[s*b+t];
        }
    }
}

#define PROTO(Ring) \
    template void Multiply \
    ( Orientation orientation, \
      Ring alpha, const SparseMatrix<Ring>& A, \
                  const Matrix<Ring>& X, \
      Ring beta,        Matrix<Ring>& Y ); \
    template void Multiply \
    ( Orientation orientation, \
      Ring alpha, const DistSparseMatrix<Ring>& A, \
                  const DistMultiVec<Ring>& X, \
      Ring beta,        DistMultiVec<Ring>& Y );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
