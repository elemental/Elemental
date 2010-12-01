      PROGRAM MAINI

      INCLUDE 'mpif.h'

      INTEGER N, NMAX, JMAX, IL, IU, TRYRAC, NZ, LDZ, OFFSET
      PARAMETER (NMAX=1000, JMAX=NMAX, LDZ=NMAX)

      DOUBLE PRECISION VL, VU, D(NMAX), E(NMAX), W(NMAX), 
     $                 Z(NMAX,JMAX), ZSUPP(2*NMAX)

      INTEGER I, J, UFILE, IERR, REQ, PROV, PID, NPROC
      PARAMETER (UFILE=1)

*     external functions
      EXTERNAL PMRRR

      REQ = MPI_THREAD_MULTIPLE
      CALL MPI_INIT_THREAD(REQ, PROV, IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, PID, IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERR)

*     Read in data
      OPEN(UFILE, FILE='Wilkinson21.data', STATUS='OLD', IOSTAT=IERR,
     $     FORM='FORMATTED')

      READ(UFILE, 10) N
      DO 100, I=1,N
         READ(UFILE, 20, END=200) D(I), E(I)
 100  CONTINUE
 200  CONTINUE

      CLOSE(UFILE)

*     Call MRRR to compute eigenpairs 3 to 18
      TRYRAC = 1
      IL = 3
      IU = 18

      CALL PMRRR('V', 'I', N, D, E, VL, VU, IL, IU, TRYRAC, 
     $           MPI_COMM_WORLD, NZ, OFFSET, W, Z, LDZ, ZSUPP, IERR)
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'Routine has failed with error', IERR
      ENDIF

*     Possibly communicate eigenvalues
*     CALL PMR_COMM_EIGVALS(MPI_COMM_WORLD, NZ, OFFSET, W)

*     Write out results; this is just a hack and will sometimes fail
*     since the output different processes may interleave
      OPEN(UFILE, FILE='result_ind.m', STATUS='UNKNOWN')
      CLOSE(UFILE)

      DO 300, I=1,NPROC
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
         IF (I .EQ. PID+1) THEN
            OPEN(UFILE, FILE='result_ind.m', STATUS='OLD', 
     $           ACCESS='APPEND')
            CALL PRTVEC(UFILE, 'W', W, NZ, OFFSET)
            CLOSE(UFILE)
         ENDIF
 300  CONTINUE
      
      DO 400, I=1,NPROC
         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
         IF (I .EQ. PID+1) THEN
            OPEN(UFILE, FILE='result_ind.m', STATUS='OLD', 
     $           ACCESS='APPEND')
            DO 500, J=1,NZ
               CALL PRTCOL(UFILE, 'Z', Z(1,J), N, J, OFFSET)
 500        CONTINUE
            CLOSE(UFILE)
         ENDIF
 400     CONTINUE

      CALL MPI_FINALIZE(IERR)

 10   FORMAT(1X, I4)
 20   FORMAT(1X, E23.17, 1X, E23.17E2)

      END




      SUBROUTINE PRTVEC(UFILE, NAME, V, N, OFFSET)

      CHARACTER NAME
      INTEGER I, N, OFFSET, UFILE
      DOUBLE PRECISION V(N)

      DO 100, I=1,N
         WRITE(UFILE,10) NAME,'(',I+OFFSET,')=',V(I),';'
 100  CONTINUE

 10   FORMAT(1X, A1, A1, I4, A2, 1X, E25.17E3, A1)
      
      END




      SUBROUTINE PRTCOL(UFILE, NAME, V, N, COL, OFFSET)

      CHARACTER NAME
      INTEGER I, N, COL, OFFSET, UFILE
      DOUBLE PRECISION V(N)

      DO 100, I=1,N
         WRITE(UFILE,10) NAME,'(',I,',',COL+OFFSET,')=',V(I),';'
 100  CONTINUE

 10   FORMAT(1X, A1, A1, I4, A1, I4, A2, 1X, E25.17E3, A1)
      
      END
