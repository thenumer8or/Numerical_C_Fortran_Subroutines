C***********************************************************************
C                                                                      *
C                 WIELANDT'S DEFLATION                                 *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SECOND MOST DOMINANT EIGENVALUE AND AN
C     ASSOCIATED EIGENVECTOR OF THE N BY N MATRIX A GIVEN AN
C     APPROXIMATION LAMBDA TO THE DOMINANT EIGENVALUE, AN
C     APPROXIMATION V TO A CORRESPONDING EIGENVECTOR AND A VECTOR X
C     BELONGING TO R**(N-1):
C
C     INPUT:   DIMENSION N; MATRIX A; APPROXIMATE EIGENVALUE LAMBDA;
C              APPROXIMATE EIGENVECTOR V BELONGING TO R**N; VECTOR X
C              BELONGING TO R**(N-1).
C
C     OUTPUT:  APPROXIMATE EIGENVALUE MU; APPROXIMATE EIGENVECTOR U OR
C              A MESSAGE THAT THE METHOD FAILS.
C
C     INITIALIZATION
      DIMENSION A(3,3),B(2,2),X(2),V(3),W(3),Y(2),VV(3)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is Wielandt Deflation.'
      WRITE(6,*) 'The array will be input from a text file in the'
      WRITE(6,*) ' order: A(1,1), A(1,2), ..., A(1,N), A(2,1),'
      WRITE(6,*) ' A(2,2), ..., A(2,N)..., A(N,1), A(N,2),'
      WRITE(6,*) ' ..., A(N,N) '
      WRITE(6,*) 'Next place the approximate eigenvector'
      WRITE(6,*) 'V(1) ... V(n), follow it with an'
      WRITE(6,*) 'approximate eigenvalue and finally with an'
      WRITE(6,*) 'initial approximate eigenvector of dimension'
      WRITE(6,*) 'n-1: X(1), ..., X(n-1)'
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
19       IF (OK) GOTO 111
         WRITE(6,*) 'Input the dimension n.'
         WRITE(6,*)
         READ(5,*) N
         IF (N .GT. 0) THEN
            M = N-1
            READ(INP,*) ((A(I,J),J=1,N),I=1,N)
            READ(INP,*) (V(I),I=1,N)
            READ(INP,*) XMU
            READ(INP,*) (X(I),I=1,M)
            OK = .TRUE.
            CLOSE(UNIT=INP)
         ELSE
            WRITE(6,*) 'The number must be a positive integer'
         ENDIF
         GOTO 19
111      OK = .FALSE.
112      IF (OK) GOTO 113
         WRITE(6,*) 'Input the tolerance for the power method.'
         WRITE(6,*) ' '
         READ(5,*) TOL
         IF (TOL .GT. 0.0) THEN
            OK = .TRUE.
         ELSE
            WRITE(6,*) 'Tolerance must be positive.'
         ENDIF
         GOTO 112
113      OK = .FALSE.
114      IF (OK) GOTO 115
         WRITE(6,*) 'Input maximum number of iterations.'
         WRITE(6,*)
         READ(5,*) NN
         IF (NN .GT. 0) THEN
            OK = .TRUE.
         ELSE
            WRITE(6,*) 'Number must be a positive integer.'
         ENDIF
         GOTO 114
      ELSE
         WRITE(6,*) 'The program will end so the input file can '
         WRITE(6,*) 'be created. '
         OK = .FALSE.
      ENDIF
115   IF(.NOT. OK) GOTO 400
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
      WRITE(OUP,*) 'WIELANDT DEFLATION'
      WRITE(OUP,5)
      WRITE(OUP,6) ((A(I,J),J=1,N),I=1,N)
      WRITE(OUP,7)
      WRITE(OUP,6) (V(I),I=1,N)
      WRITE(OUP,8) XMU
      WRITE(OUP,9) (X(I),I=1,M)
C     STEP 1
      I=1
      AMAX=ABS(V(1))
      DO 100 J=2,N
           IF(ABS(V(J)).GT.AMAX) THEN
                I=J
                AMAX=ABS(V(J))
           END IF
100   CONTINUE
      IF (AMAX.LT.1.0E-20) THEN
           WRITE(OUP,10)
           GOTO 400
      END IF
      I1=I-1
      N1=N-1
C     STEP 2
      IF (I .NE. 1) THEN
           DO 20 K=1,I1
                DO 20 J=1,I1
20         B(K,J) = A(K,J)-V(K)*A(I,J)/V(I)
      END IF
C     STEP 3
      IF(I.NE.1 .AND. I.NE.N) THEN
           DO 30 K=I,N1
           DO 30 J=1,I1
                B(K,J) = A(K+1,J)-V(K+1)*A(I,J)/V(I)
30         B(J,K) = A(J,K+1)-V(J)*A(I,K+1)/V(I)
      END IF
C     STEP 4
      IF (I.NE.N) THEN
           DO 40 K=I,N1
                DO 40 J=I,N1
40         B(K,J) = A(K+1,J+1)-V(K+1)*A(I,J+1)/V(I)
      END IF
      WRITE(OUP,11)
      WRITE(OUP,12) ((B(L1,L2),L2=1,M),L1=1,M)
C     STEP 5
      CALL POWER(M,B,X,Y,YMU,NN,TOL)
C     Y IS USED IN PLACE OF W'
C     STEP 6
      IF (I.NE.1) THEN
           DO 50 K=1,I1
50         W(K) = Y(K)
      END IF
C     STEP 7
      W(I) = 0
C     STEP 8
      IF (I.NE.N) THEN
           I2 = I+1
           DO 60 K=I2,N
60         W(K) = Y(K-1)
      END IF
C     STEP 9
      DO 70 K=1,N
           S = 0
           DO 80 J=1,N
80         S = S+A(I,J)*W(J)
           S = S/V(I)
C     COMPUTE EIGENVECTOR
C     VV IS USED IN PLACE OF U HERE
70    VV(K) = (YMU-XMU)*W(K)+S*V(K)
      WRITE(OUP,13) YMU
      WRITE(OUP,14) (Y(I),I=1,M)
      WRITE(OUP,15) (VV(I),I=1,N)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
5     FORMAT(1X,'ORIGINAL MATRIX A')
6     FORMAT(1X,3(3X,F10.5))
7     FORMAT(1X,'APPROXIMATE EIGENVECTOR')
8     FORMAT(1X,'APPROX. DOMINANT EIGENVALUE',3X,E15.8)
9     FORMAT(1X,'INITIAL APPROX.',2(3X,F10.5))
10    FORMAT(1X,'INPUT ERROR')
11    FORMAT(1X,'REDUCED MATRIX B')
12    FORMAT(1X,2(3X,E15.8))
13    FORMAT(1X,'APPROX. TO SECOND EIGENVALUE',3X,E15.8)
14    FORMAT(1X,'EIGENVECTOR FOR B',2(3X,E15.8))
15    FORMAT(1X,'EIGENVECTOR FOR A',3(3X,E15.8))
      END
C
      SUBROUTINE POWER(N,A,X,Y,YMU,NN,TOL)
C     SEE ALGORITHM 9.1
C     NOTATION: N X N MATRIX A, ITERATION VECTOR Y, NORMALIZED
C     APPROX. EIGENVECTOR X, APPROX. EIGENVALUE XMU
      DIMENSION A(N,N),X(N),Y(N)
      K = 1
      LP = 1
      AMAX = ABS(X(1))
      DO 10 I=2,N
           IF(ABS(X(I)).GT.AMAX) THEN
                AMAX = ABS(X(I))
                LP = I
           END IF
10    CONTINUE
      DO 20 I=1,N
20    X(I)=X(I)/AMAX
110   IF (K.GT.NN) GOTO 200
           DO 30 I=1,N
                Y(I) = 0.0
                DO 30 J=1,N
30         Y(I) = Y(I)+A(I,J)*X(J)
           YMU=Y(LP)
           LP = 1
           AMAX = ABS(Y(1))
           DO 40 I=2,N
                IF(ABS(Y(I)).GT.AMAX) THEN
                     AMAX = ABS(Y(I))
                     LP = I
                END IF
40        CONTINUE
           IF (AMAX.LE.TOL) THEN
                WRITE(6,100)
                STOP
           ENDIF
           ERR=0.0
           DO 50 I=1,N
              T=Y(I)/Y(LP)
              IF(ABS(X(I)-T).GT.ERR) ERR=ABS(X(I)-T)
50         X(I) = T
           IF(ERR.LE.TOL) THEN
                   DO 60 I=1,N
60                 Y(I)=X(I)
                   RETURN
           END IF
           K=K+1
      GOTO 110
200   WRITE(6,101)
      STOP
100   FORMAT(1X,'ZERO VECTOR-ERROR')
101   FORMAT(1X,'DIVERGENCE-ERROR')
      END
