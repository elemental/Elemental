/*
   Copyright (c) 2009-2012, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
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

   Authors:
   This interface is mainly due to Martin Schatz, but it was put into its
   current form by Jack Poulson.
*/

namespace elem {

namespace axpy_type_wrapper {
enum AxpyType { LOCAL_TO_GLOBAL, GLOBAL_TO_LOCAL };
}
using namespace axpy_type_wrapper;

template<typename T,typename Int=int>
class AxpyInterface
{   
public:
    AxpyInterface();
    ~AxpyInterface();

    AxpyInterface( AxpyType type,       DistMatrix<T,MC,MR>& Z );
    AxpyInterface( AxpyType type, const DistMatrix<T,MC,MR>& Z ); 

    void Attach( AxpyType type,       DistMatrix<T,MC,MR>& Z ); 
    void Attach( AxpyType type, const DistMatrix<T,MC,MR>& Z ); 

    void Axpy( T alpha,       Matrix<T>& Z, Int i, Int j );
    void Axpy( T alpha, const Matrix<T>& Z, Int i, Int j );

    void Detach();

private:
    static const Int 
        DATA_TAG        =1, 
        EOM_TAG         =2, 
        DATA_REQUEST_TAG=3, 
        DATA_REPLY_TAG  =4;

    bool attachedForLocalToGlobal_, attachedForGlobalToLocal_;
    byte sendDummy_, recvDummy_;
    DistMatrix<T,MC,MR>* localToGlobalMat_;
    const DistMatrix<T,MC,MR>* globalToLocalMat_;

    std::vector<bool> sentEomTo_, haveEomFrom_;
    std::vector<byte> recvVector_;
    std::vector<mpi::Request> eomSendRequests_;

    std::vector<std::deque<std::vector<byte> > >
        dataVectors_, requestVectors_, replyVectors_;
    std::vector<std::deque<bool> > 
        sendingData_, sendingRequest_, sendingReply_;
    std::vector<std::deque<mpi::Request> > 
        dataSendRequests_, requestSendRequests_, replySendRequests_;

    // Check if we are done with this attachment's work
    bool Finished();

    // Progress functions
    void UpdateRequestStatuses();
    void HandleEoms();
    void HandleLocalToGlobalData();
    void HandleGlobalToLocalRequest();
    void StartSendingEoms();
    void FinishSendingEoms();

    void AxpyLocalToGlobal( T alpha, const Matrix<T>& X, Int i, Int j );
    void AxpyGlobalToLocal( T alpha,       Matrix<T>& Y, Int i, Int j );

    Int ReadyForSend
    ( Int sendSize,
      std::deque<std::vector<byte> >& sendVectors,
      std::deque<mpi::Request>& requests, 
      std::deque<bool>& requestStatuses );
};

} // namespace elem
