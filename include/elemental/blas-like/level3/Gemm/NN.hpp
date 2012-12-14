/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Normal Normal Gemm that avoids communicating the matrix A.
template<typename T>
inline void
GemmNNA
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNNA");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNNA: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);
    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);

    // Temporary distributions
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);
    DistMatrix<T,MC,STAR> D1_MC_STAR(g);

    B1_VR_STAR.AlignWith( A );
    B1Trans_STAR_MR.AlignWith( A );
    D1_MC_STAR.AlignWith( A );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
    while( BR.Width() > 0 )
    {
        LockedRepartitionRight
        ( BL, /**/     BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/     CR,
          C0, /**/ C1, C2 );

        Zeros( C1.Height(), C1.Width(), D1_MC_STAR );
        //--------------------------------------------------------------------//
        B1_VR_STAR = B1;
        B1Trans_STAR_MR.TransposeFrom( B1_VR_STAR );

        // D1[MC,*] := alpha A[MC,MR] B1[MR,*]
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, A, B1Trans_STAR_MR, T(0), D1_MC_STAR );

        // C1[MC,MR] += scattered result of D1[MC,*] summed over grid rows
        C1.SumScatterUpdate( T(1), D1_MC_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm that avoids communicating the matrix B.
template<typename T>
inline void 
GemmNNB
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNNB");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNNB: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);

    // Temporary distributions
    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( B );
    D1Trans_MR_STAR.AlignWith( B );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    PartitionDown
    ( C, CT,
         CB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        RepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );

        Zeros( C1.Width(), C1.Height(), D1Trans_MR_STAR );
        //--------------------------------------------------------------------//
        A1_STAR_MC = A1; // A1[*,MC] <- A1[MC,MR]

        // D1^T[MR,* ] := alpha B^T[MR,MC] A1^T[MC,* ]
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, alpha, B, A1_STAR_MC, T(0), D1Trans_MR_STAR );

        C1.TransposeSumScatterUpdate( T(1), D1Trans_MR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
 
        SlidePartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}                     

// Normal Normal Gemm that avoids communicating the matrix C.
template<typename T>
inline void 
GemmNNC
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNNC");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNNC: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);         
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    // Temporary distributions
    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g); 

    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    // Start the algorithm
    Scale( beta, C );
    LockedPartitionRight( A, AL, AR, 0 ); 
    LockedPartitionDown
    ( B, BT, 
         BB, 0 ); 
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight( AL, /**/ AR,
                                A0, /**/ A1, A2 );

        LockedRepartitionDown( BT,  B0,
                              /**/ /**/
                                    B1, 
                               BB,  B2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1; 
        B1Trans_MR_STAR.TransposeFrom( B1 );

        // C[MC,MR] += alpha A1[MC,*] (B1^T[MR,*])^T
        //           = alpha A1[MC,*] B1[*,MR]
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, A1_MC_STAR, B1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight( AL,     /**/ AR,
                                   A0, A1, /**/ A2 );

        SlideLockedPartitionDown( BT,  B0,
                                       B1,
                                 /**/ /**/
                                  BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// Normal Normal Gemm for panel-panel dot products. 
template<typename T>
inline void 
GemmNNDot
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNNDot");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() ||
        B.Width()  != C.Width()  ||
        A.Width()  != B.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal GemmNNDot: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    if( A.Height() > B.Width() )
    {
        // Matrix views
        DistMatrix<T> AT(g), AB(g),
                      A0(g), A1(g), A2(g);         
        DistMatrix<T> BL(g),  B0(g),
                      BR(g),  B1(g),
                              B2(g);
        DistMatrix<T> CT(g), C0(g), C1L(g), C1R(g),
                      CB(g), C1(g), C10(g), C11(g), C12(g),
                             C2(g);

        // Temporary distributions
        DistMatrix<T,STAR,VC> A1_STAR_VC(g);
        DistMatrix<T,VC,STAR> B1_VC_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        // Star the algorithm
        Scale( beta, C );
        LockedPartitionDown
        ( A, AT,
             AB, 0 );
        PartitionDown
        ( C, CT,
             CB, 0 );
        while( AB.Height() > 0 )
        {
            LockedRepartitionDown
            ( AT,  A0,
             /**/ /**/
                   A1,
              AB,  A2 );

            RepartitionDown
            ( CT,  C0,
             /**/ /**/
                   C1,
              CB,  C2 );

            A1_STAR_VC = A1; 
            B1_VC_STAR.AlignWith( A1_STAR_VC );

            LockedPartitionRight( B, BL, BR, 0 );
            PartitionRight( C1, C1L, C1R, 0 );
            while( BR.Width() > 0 )
            {
                LockedRepartitionRight
                ( BL, /**/ BR,
                  B0, /**/ B1, B2 );

                RepartitionRight
                ( C1L, /**/ C1R,
                  C10, /**/ C11, C12 );

                Zeros( C11.Height(), C11.Width(), C11_STAR_STAR );
                //------------------------------------------------------------//
                B1_VC_STAR = B1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VC, B1_VC_STAR, T(0), C11_STAR_STAR );
                C11.SumScatterUpdate( T(1), C11_STAR_STAR );
                //------------------------------------------------------------//

                SlideLockedPartitionRight
                ( BL,     /**/ BR,
                  B0, B1, /**/ B2 );

                SlidePartitionRight
                ( C1L,      /**/ C1R,
                  C10, C11, /**/ C12 );
            }
            B1_VC_STAR.FreeAlignments();

            SlideLockedPartitionDown
            ( AT,  A0,
                   A1,
             /**/ /**/
              AB,  A2 );

            SlidePartitionDown
            ( CT,  C0,
                   C1,
             /**/ /**/
              CB,  C2 );
        }
    }
    else
    {
        // Matrix views
        DistMatrix<T> AT(g), AB(g),
                      A0(g), A1(g), A2(g);         
        DistMatrix<T> BL(g),  B0(g),
                      BR(g),  B1(g),
                              B2(g);
        DistMatrix<T> 
            CL(g), CR(g),         C1T(g),  C01(g),
            C0(g), C1(g), C2(g),  C1B(g),  C11(g),
                                           C21(g);

        // Temporary distributions
        DistMatrix<T,STAR,VR> A1_STAR_VR(g);
        DistMatrix<T,VR,STAR> B1_VR_STAR(g);
        DistMatrix<T,STAR,STAR> C11_STAR_STAR(g);

        // Star the algorithm
        Scale( beta, C );
        LockedPartitionRight( B, BL, BR, 0 );
        PartitionRight( C, CL, CR, 0 );
        while( BR.Width() > 0 )
        {
            LockedRepartitionRight
            ( BL, /**/ BR,
              B0, /**/ B1, B2 );

            RepartitionRight
            ( CL, /**/ CR,
              C0, /**/ C1, C2 );

            B1_VR_STAR = B1;
            A1_STAR_VR.AlignWith( B1_VR_STAR );

            LockedPartitionDown
            ( A, AT,
                 AB, 0 );
            PartitionDown
            ( C1, C1T,
                  C1B, 0 );
            while( AB.Height() > 0 )
            {
                LockedRepartitionDown
                ( AT,  A0,
                 /**/ /**/
                       A1,
                  AB,  A2 );

                RepartitionDown
                ( C1T,  C01,
                 /***/ /***/
                        C11,
                  C1B,  C21 );

                Zeros( C11.Height(), C11.Width(), C11_STAR_STAR );
                //------------------------------------------------------------//
                A1_STAR_VR = A1;
                LocalGemm
                ( NORMAL, NORMAL, 
                  alpha, A1_STAR_VR, B1_VR_STAR, T(0), C11_STAR_STAR );
                C11.SumScatterUpdate( T(1), C11_STAR_STAR );
                //------------------------------------------------------------//

                SlideLockedPartitionDown
                ( AT,  A0,
                       A1,
                 /**/ /**/
                  AB,  A2 );

                SlidePartitionDown
                ( C1T,  C01,
                        C11,
                 /***/ /***/
                  C1B,  C21 );
            }
            A1_STAR_VR.FreeAlignments();

            SlideLockedPartitionRight
            ( BL,     /**/ BR,
              B0, B1, /**/ B2 ); 

            SlidePartitionRight
            ( CL,     /**/ CR,
              C0, C1, /**/ C2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
GemmNN
( T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::GemmNN");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
#endif
    const int m = C.Height();
    const int n = C.Width();
    const int k = A.Width();
    const float weightTowardsC = 2.0;
    const float weightAwayFromDot = 10.0;

    if( weightAwayFromDot*m <= k && weightAwayFromDot*n <= k )
    {
        GemmNNDot( alpha, A, B, beta, C );
    }
    else if( m <= n && weightTowardsC*m <= k )
    {
        GemmNNB( alpha, A, B, beta, C );    
    }
    else if( n <= m && weightTowardsC*n <= k )
    {
        GemmNNA( alpha, A, B, beta, C );
    }
    else
    {
        GemmNNC( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
