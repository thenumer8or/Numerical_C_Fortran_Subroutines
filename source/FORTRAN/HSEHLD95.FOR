C***********************************************************************
C                                                                      *
C                    HOUSEHOLDER'S METHOD                              *
C                                                                      *
C***********************************************************************
C
C
C
C     TO OBTAIN A SYMMETRIC TRIDIAGONAL MATRIX A(N-1) SIMILAR
C     TO THE SYMMETRIC MATRIX A=A(1), CONSTRUCT THE FOLLOWING
C     MATRICES A(2),A(3),...,A(N-1) WHERE A(K) = A(I,J)**K, FOR
C     EACH K = 1,2,...,N-1:
C
C     INPUT:   DIMENSION N; MATRIX A.
C
C     OUTPUT:  A(N-1) (COULD OVER-WRITE A).
C
C     INITIALIZATION
      DIMENSION A(10,10),U(10),V(10),Z(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Householder Method.'
      WRITE(6,*) 'The symmetric array A will be input from a text'
      WRITE(6,*) 'in the order:'
      WRITE(6,*) 'A(1,1), A(1,2), A(1,3), ..., A(1,n)'
      WRITE(6,*) '        A(2,2), A(2,3), ..., A(2,n)'
      WRITE(6,*) '                A(3,3), ..., A(3,n)'
      WRITE(6,*) '                        ..., A(n,n)'
      WRITE(6,*) 'Place as many entries as desired on each line,'
      WRITE(6,*) ' but separate entries with at least one blank.'
      OK = .FALSE.
      WRITE(6,*) 'Has the input file been created?'
      WRITE(6,*) 'Enter Y or N - letter within quotes '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF (( AA .EQ. 'Y' ) .OR.( AA .EQ. 'y' )) THEN
         WRITE(6,*) 'Input the file name in the form - '
         WRITE(6,*) 'drive:name.ext  contained in quotes'
         WRITE(6,*) 'as example:   ''A:DATA.DTA'' '
         WRITE(6,*) ' '
         READ(5,*)  NAME
         INP = 4
         OPEN(UNIT=INP,FILE=NAME,ACCESS='SEQUENTIAL')
         OK = .FALSE.
9        IF (OK) GOTO 11
         WRITE(6,*) 'Input the dimension n.'
         WRITE(6,*)
         READ(5,*) N
         IF (N .GT. 0) THEN
            READ(INP,*) ((A(I,J),J=I,N),I=1,N)
            DO 41 I=1,N
               DO 41 J=I,N
                  A(J,I) = A(I,J)
41          CONTINUE
            OK = .TRUE.
            CLOSE(UNIT=INP)
         ELSE
            WRITE(6,*) 'The number must be a positive integer'
         ENDIF
         GOTO 9
      ELSE
         WRITE(6,*) 'The program will end so the input file can '
         WRITE(6,*) 'be created. '
         OK = .FALSE.
      ENDIF
11    IF(.NOT. OK) GOTO 400
      WRITE(6,*) 'Select output destinations: '
      WRITE(6,*) '1. Screen '
      WRITE(6,*) '2. Text file '
      WRITE(6,*) 'Enter 1 or 2 '
      WRITE(6,*) ' '
      READ(5,*) FLAG
      IF ( FLAG .EQ. 2 ) THEN
         WRITE(6,*) 'Input the file name in the form - '
         WRITE(6,*) 'drive:name.ext'
         WRITE(6,*) 'with the name contained within quotes'
         WRITE(6,*) 'as example:   ''A:OUTPUT.DTA'' '
         WRITE(6,*) ' '
         READ(5,*) NAME1
         OUP = 3
         OPEN(UNIT=OUP,FILE=NAME1,STATUS='NEW')
      ELSE
         OUP = 6
      ENDIF
      WRITE(OUP,*) 'HOUSEHOLDER METHOD'
      WRITE(OUP,3)
      WRITE(OUP,4) ((A(I,J),J=1,N),I=1,N)
C     STEP 1
      N2 = N-2
      N1 = N-1
      DO 10 K=1,N2
           I2 = K+1
           I3 = K+2
C          STEP 2
           Q = 0
           DO 20 J=I2,N
20              Q = Q+A(J,K)**2
           Q = SQRT(Q)
C          STEP 3
           IF(ABS(A(K+1,K)).LE.1.0E-40) THEN
                S=Q
           ELSE
                S = Q*(A(K+1,K)/ABS(A(K+1,K)))
           END IF
C          STEP 4
           RSQ = S*S+S*A(K+1,K)
C          NOTE: RSQ = 2*R**2
C          STEP 5
           V(K) = 0
C          NOTE: V(1)=...=V(K-1)=0, BUT ARE NOT NEEDED
           V(K+1) = A(K+1,K)+S
           IF(K.LT.N-1)THEN
                DO 30 J=I3,N
30              V(J) = A(J,K)
           END IF
C          NOTE: W=V/SQRT(2*RSQ)=V/(2*R)
C          STEP 6
           DO 40 J=K,N
                U(J) = 0
                DO 50 L=I2,N
50                   U(J) = U(J)+A(J,L)*V(L)
40              U(J) = U(J)/RSQ
C          NOTE: U=A(K)*V/RSQ=A(K)*V/(2*R**2)
C          STEP 7
           PROD = 0
           DO 60 J=I2,N
60             PROD = PROD+V(J)*U(J)
C          NOTE: PROD = V ** T *U = V**T *A(K)*V/(2*R**2)
C          STEP 8
           QUO = PROD/(2*RSQ)
           DO 70 J=K,N
70              Z(J) = U(J)-QUO*V(J)
C          NOTE: Z=U-V**TUV/(2RSQ)=U-V**TUV/(4R**2)=U-WW**TU
C                 =A(K)*W/R-WW**T(A(K)W/R)
C          STEP 9
           DO 80 L=I2,N1
                L1 = L+1
C               STEP 10
                DO 90 J=L1,N
                     A(L,J) = A(L,J)-V(L)*Z(J)-V(J)*Z(L)
90                   A(J,L) = A(L,J)
C               STEP 11
80              A(L,L) = A(L,L)-2*V(L)*Z(L)
C          STEP 12
           A(N,N) = A(N,N)-2*V(N)*Z(N)
           IF (K.LT.N-1) THEN
C               STEP 13
                DO 100 J=I3,N
                     A(K,J) = 0
100                  A(J,K) = 0
           END IF
C          STEP 14
           A(K+1,K) = A(K+1,K)-V(K+1)*Z(K)
           A(K,K+1) = A(K+1,K)
           KK = K+1
           WRITE(OUP,5) KK
           WRITE(OUP,4)((A(I,J),J=1,N),I=1,N)
C     NOTE: THE OTHER ELEMENTS OF A(K+1) ARE THE SAME AS A(K)
C     A(K+1)=A(K)-VZ**T-ZV**T=(I-2WW**T)A(K)(I-2WW**T)
10    CONTINUE
C     STEP 15
C     OUTPUT HAS ALREADY BEEN PERFORMED
C     THE PROCESS IS COMPLETE, A(N-1) IS SYMMETRIC, TRIDIAGONAL AND
C     SIMILAR TO A
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'ORIGINAL MATRIX A=A(1)')
4     FORMAT(1X,4(3X,E15.8))
5     FORMAT(1X,'THE MATRIX A(',I2,') EQUALS')
      END
