C***********************************************************************
C                                                                      *
C                           QR METHOD                                  *
C                                                                      *
C***********************************************************************
C
C
C
C     TO OBTAIN THE EIGENVALUES OF A SYMMETRIC, TRIDIAGONAL N BY N
C     MATRIX:
C
C                 A(1)  B(2)
C                 B(2)  A(2)  B(3)
C                      .     .      .
C                         .      .     .
C                           .      .     .
C                         B(N-1)  A(N)  B(N)
C                                 B(N)  A(N)
C
C     INPUT:   N; A(1),...,A(N) (DIAGONAL OF A); B(2),...,B(N)
C              (OFF-DIAGONAL OF A); MAXIMUM NUMBER OF ITERATIONS M;
C              TOLERANCE TOL
C
C     OUTPUT:  EIGENVALUES OF A OR RECOMMENDED SPLITTING OF A
C              OR A MESSAGE THAT THE MAXIMUM NUMBER OF
C              ITERATIONS WERE EXCEEDED.
C
      DIMENSION A(9),B(9),C(9),D(9),Q(9),R(9),S(9),X(9),Y(9),Z(9)
C     INPUT N, A(1),...,A(N), B(2),...,B(N), MM, TOL
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the QR Method.'
      WRITE(6,*) 'The tridiagonal symmetric array A will be input'
      WRITE(6,*) 'from a text file in the order:'
      WRITE(6,*) '(diagonal):     A(1), A(2), ..., A(n)'
      WRITE(6,*) '(subdiagonal):  B(2), B(3), ..., B(n)'
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
            READ(INP,*) (A(I),I=1,N), (B(I),I=2,N)
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
         READ(5,*) MM
         IF (MM .GT. 0) THEN
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
      WRITE(OUP,*) 'QR METHOD'
      ZERO = 1.0E-30
C     STEP 1
C     SET THE ACCUMULATED SHIFT TO ZERO
      SHIFT = 0
      K=1
C     STEP 2
100   IF (K.GT.MM) GOTO 200
           WRITE(OUP,3) K
           WRITE(OUP,4) (A(I),I=1,N), (B(I),I=2,N)
C          TEST FOR CONVERGENCE / DEFLATION
C          TEST FOR POSSIBLE SPLIT
C          STEP 3
           IF (ABS(B(N)).LE.TOL) THEN
                A(N)=A(N)+SHIFT
                WRITE(OUP,6) A(N)
                N=N-1
           END IF
           IF(ABS(B(2)).LE.TOL) THEN
                A(1)=A(1)+SHIFT
                WRITE(OUP,6) A(1)
                N=N-1
                A(1)=A(2)
                DO 10 I=2,N
                     A(I)=A(I+1)
10              B(I)=B(I+1)
           END IF
           M=N-1
           IF(M.GE.2) THEN
                DO 20 I=3,M
                     IF(ABS(B(I)).LE.TOL) THEN
                          WRITE(OUP,5)
                          GOTO 400
                     END IF
20              CONTINUE
           END IF
C          STEP 4
C          COMPUTE SHIFT
           B1=-(A(N-1)+A(N))
           C1=A(N)*A(N-1)-B(N)*B(N)
           IF (ABS(C1).LE.ZERO) THEN
                WRITE(OUP,7)
C               GOTO 400
           END IF
           D1=B1*B1-4*C1
           IF(D1.LT.0.0) THEN
                WRITE(OUP,8)
                GOTO 400
           END IF
           D1 = SQRT(D1)
           IF(B1.GT.0.0) THEN
                X1=-2*C1/(B1+D1)
                X2=-(B1+D1)/2
           ELSE
                X1=(D1-B1)/2
                X2=2*C1/(D1-B1)
           END IF
C          IF N IS 2 THEN WE HAVE COMPUTED THE 2 EIGENVALUES
           IF(N.EQ.2) THEN
                X1=X1+SHIFT
                X2=X2+SHIFT
                WRITE(OUP,6) X1
                WRITE(OUP,6) X2
                GOTO 400
           END IF
C          STEP 5
C          ACCUMULATE SHIFT
           IF(ABS(A(N)-X1).GE.ABS(A(N)-X2)) X1=X2
           SHIFT=SHIFT+X1
           WRITE(OUP,9) X1
C          STEP 6
C          PERFORM SHIFT
           DO 30 I=1,N
30         D(I)=A(I)-X1
C          STEP 7
C          COMPUTE R(K)
           X(1)=D(1)
           Y(1)=B(2)
           DO 40 J=2,N
                Z(J-1)=SQRT(X(J-1)*X(J-1)+B(J)*B(J))
                C(J)=X(J-1)/Z(J-1)
                S(J)=B(J)/Z(J-1)
                Q(J-1)=C(J)*Y(J-1)+S(J)*D(J)
                X(J)=C(J)*D(J)-S(J)*Y(J-1)
                IF(J.NE.N) THEN
                     R(J-1)=S(J)*B(J+1)
                     Y(J)=C(J)*B(J+1)
                END IF
40         CONTINUE
           M1=N-2
C          STEP 8
C          COMPUTE NEW A
           Z(N)=X(N)
           A(1)=S(2)*Q(1)+C(2)*Z(1)
           B(2)=S(2)*Z(2)
           IF(N.GT.2) THEN
                DO 50 J=2,M
                     A(J)=S(J+1)*Q(J)+C(J+1)*C(J)*Z(J)
50              B(J+1)=S(J+1)*Z(J+1)
           END IF
           A(N)=C(N)*X(N)
           K=K+1
      GOTO 100
200   WRITE(OUP,12)
C     STEP 9
C     THE PROCESS IF COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,I2)
4     FORMAT(1X,'Vectors A-B',7(2X,E15.8))
5     FORMAT(1X,'TIME TO SPLIT')
6     FORMAT(1X,'EIGENVALUE = ',E15.8,' DEFLATE AND MOVE ON')
7     FORMAT(1X,'NEARLY SINGULAR')
8     FORMAT(1X,'COMPLEX ROOTS')
9     FORMAT(1X,'SHIFT = ',E15.8)
11    FORMAT(7(2X,E15.8))
12    FORMAT(1X,'FAILURE')
      END
