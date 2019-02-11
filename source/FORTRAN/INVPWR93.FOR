C***********************************************************************
C                                                                      *
C               INVERSE POWER METHOD                                   *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE AN EIGENVALUE AND AN ASSOCIATED EIGENVECTOR OF THE
C     N BY N MATRIX A GIVEN A NONZERO VECTOR X:
C
C     INPUT:   DIMENSION N; MATRIX A; VECTOR X; TOLERANCE TOL;
C              MAXIMUM NUMBER OF ITERATIONS N.
C
C     OUTPUT:  APPROXIMATE EIGENVALUE MU; APPROXIMATE EIGENVECTOR
C              X OR A MESSAGE THAT THE MAXIMUM NUMBER OF ITERATIONS
C              WAS EXCEEDED.
C
      DIMENSION A(3,3),X(3),Y(3),NROW(3),B(3)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Inverse Power Method.'
      WRITE(6,*) 'The array will be input from a text file in the'
      WRITE(6,*) ' order: A(1,1), A(1,2), ..., A(1,N), A(2,1),'
      WRITE(6,*) ' A(2,2), ..., A(2,N)..., A(N,1), A(N,2),'
      WRITE(6,*) ' ..., A(N,N) '
      WRITE(6,*) 'Place as many entries as desired on each line,'
      WRITE(6,*) ' but separate entries with at least one blank.'
      WRITE(6,*) 'The initial approximation should follow in the'
      WRITE(6,*) 'same format.'
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
            READ(INP,*) ((A(I,J), J=1,N),I=1,N)
            READ(INP,*) (X(I),I=1,N)
            OK = .TRUE.
            CLOSE(UNIT=INP)
         ELSE
            WRITE(6,*) 'The number must be a positive integer'
         ENDIF
         GOTO 19
111      OK = .FALSE.
112      IF (OK) GOTO 13
         WRITE(6,*) 'Input the tolerance.'
         WRITE(6,*) ' '
         READ(5,*) TOL
         IF (TOL .GT. 0.0) THEN
            OK = .TRUE.
         ELSE
            WRITE(6,*) 'Tolerance must be positive.'
         ENDIF
         GOTO 112
13       OK = .FALSE.
14       IF (OK) GOTO 15
         WRITE(6,*) 'Input maximum number of iterations.'
         WRITE(6,*)
         READ(5,*) NN
         IF (NN .GT. 0) THEN
            OK = .TRUE.
         ELSE
            WRITE(6,*) 'Number must be a positive integer.'
         ENDIF
         GOTO 14
      ELSE
         WRITE(6,*) 'The program will end so the input file can '
         WRITE(6,*) 'be created. '
         OK = .FALSE.
      ENDIF
15    IF(.NOT. OK) GOTO 400
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
      WRITE(OUP,*) 'INVERSE POWER METHOD'
C     STEP 1
C     Q COULD BE INPUT INSTEAD OF COMPUTED BY DELETING THE NEXT 7 STEPS
      Q=0
      S = 0
      DO 50 I=1,N
           S= S + X(I)*X(I)
           DO 50 J=1,N
50         Q = Q + A(I,J)*X(I)*X(J)
      Q = Q/S
C     STEP 2
      K =1
      WRITE(OUP,4)
      WRITE(OUP,5) ((A(I,J),J=1,N),I=1,N)
      WRITE(OUP,*) 'INPUT VECTOR'
      J=0
      WRITE(OUP,6) J,(X(I),I=1,N)
      WRITE(OUP,7) Q
C     FORM MATRIX A - Q*I
      DO 10 I=1,N
10         A(I,I) = A(I,I)-Q
C     CALL SUBROUTINE TO COMPUTE MULTIPLIERS M(I,J) AND UPPER TRIANGULAR
C     MATRIX FOR MATRIX A USING GAUSS ELIMINATION WITH MAXIMAL COLUMN
C     PIVOTING--NROW HOLDS THE ORDERING OF ROWS FOR INTERCHANGES
      CALL MULTIP(N,A,NROW)
C     STEP 3
      LP = 1
      DO 20 I=2,N
20         IF(ABS(X(I)).GT.ABS(X(LP))) LP =I
C     STEP 4
      AMAX = X(LP)
      DO 30 I=1,N
30         X(I) = X(I)/AMAX
C     STEP 5
100   IF (K.GT.NN) GOTO 200
C          STEP 6
           DO 35 I=1,N
35         B(I)=X(I)
C          SUBROUTINE SOLVE RETURNS THE SOLUTION OF (A-Q*I)Y=B IN Y
           CALL SOLVE (N,A,NROW,B,Y)
C          STEP 7
           YMU = Y(LP)
C          STEP 8 AND 9
           LP=1
           DO 40 I=2,N
40              IF(ABS(Y(I)).GT.ABS(Y(LP))) LP=I
           AMAX = Y(LP)
           ERR=0.0
           DO 60 I=1,N
                T=Y(I)/AMAX
                IF(ABS(X(I)-T).GT.ERR) ERR=ABS(X(I)-T)
60              X(I)=T
           WRITE(OUP,*) 'ITER. NO. AND VECTOR'
           WRITE(OUP,6) K,(X(I),I=1,N)
           WRITE(OUP,8) YMU,ERR
C          STEP 10
           IF (ERR.LT.TOL) THEN
                YMU=1/YMU+Q
                WRITE(OUP,3) YMU
C               PROCEDURE COMPLETED SUCCESSFULLY
                GOTO 400
           END IF
C          STEP 11
           K = K+1
      GOTO 100
C     STEP 12
200   CONTINUE
      WRITE (OUP,9)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'APPROXIMATE EIGENVALUE IS ',E15.8)
4     FORMAT(3X,' MATRIX A ')
5     FORMAT(3(1X,E15.8))
6     FORMAT(1X,I2,3(1X,E15.8))
8     FORMAT(1X,'APPROX. EIGENVALUE',1X,E15.8,' ERROR',E15.8)
7     FORMAT(3X,'Q = ',E15.8)
9     FORMAT(1X,'NO COVERGENCE')
      END
C
      SUBROUTINE MULTIP (N,A,NROW)
C
C     SUBROUTINE MULTIP COMPUTES THE MULTIPLIERS AND ROW ORDERING FOR
C     GAUSS ELIMINATION WITH MAXIMAL COLUMN PIVOTING-THE MULTIPLIERS
C     AND UPPER TRIANGULAR FORM ARE RETURNED IN A AND THE ROW ORDERING
C     IS RETURNED IN NROW
C
      DIMENSION A(N,N),NROW(N)
      DO 1 I=1,N
1     NROW(I) = I
      M = N-1
      DO 2 I=1,M
           IMAX  = I
           J = I+1
           DO 3 IP=J,N
           L1 = NROW(IMAX)
           L2 = NROW(IP)
3          IF(ABS(A(L2,I)).GT.ABS(A(L1,I))) IMAX = IP
           IF(ABS(A(NROW(IMAX),I)).LE.1.0E-30) THEN
                WRITE(6,100)
                STOP
           END IF
           JJ = NROW(I)
           NROW(I) = NROW(IMAX)
           NROW(IMAX) = JJ
           I1 = NROW(I)
           DO 4 JJ=J,N
                J1 = NROW(JJ)
                A(J1,I) = A(J1,I)/A(I1,I)
                DO 5 K=J,N
5               A(J1,K) = A(J1,K)-A(J1,I)*A(I1,K)
4          CONTINUE
2     CONTINUE
      RETURN
100   FORMAT(1X,'MULTIP FAILS')
      END
C
      SUBROUTINE SOLVE(N,A,NROW,X,Y)
C
C     SUBROUTINE SOLVE ACCEPTS THE VECTOR X, THE MULTIPLIERS AND UPPER
C     TRIANGULAR FORM FOR A, THE ROW ORDERING NROW, SOLVES A*Y = X, AND
C     RETURNS THE SOLUTION IN Y
C
      DIMENSION A(N,N),NROW(N),X(N),Y(N)
      M = N-1
      DO 2 I=1,M
           J=I+1
           I1 = NROW(I)
           DO 4 JJ=J,N
                J1 = NROW(JJ)
4          X(J1) = X(J1)-A(J1,I)*X(I1)
2     CONTINUE
      IF(ABS(A(NROW(N),N)).LE.1.0E-30) THEN
           WRITE(6,100)
           STOP
      END IF
      N1 = NROW(N)
      Y(N) = X(N1)/A(N1,N)
      L = N-1
      DO 15 K=1,L
           J = L-K+1
           JJ = J+1
           N2 = NROW(J)
           Y(J)=X(N2)
           DO 16 KK=JJ,N
16         Y(J) = Y(J)-A(N2,KK)*Y(KK)
15    Y(J) = Y(J)/A(N2,J)
      RETURN
100   FORMAT(1X,'SOLVE FAILS')
      END
