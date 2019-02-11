C***********************************************************************
C                                                                      *
C                  STEEPEST DESCENT METHOD                             *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE A SOLUTION P TO THE MINIMIZATION PROBLEM
C                 G(P) = MIN ( G(X) : X IN N-DIM )
C     GIVEN AN INITIAL APPROXIMATION X:
C
C     INPUT NUMBER N OF VARIABLES; INITIAL APPROXIMATION X;
C           TOLERANCE TOL; MAXIMUM NUMBER OF ITERATIONS N.
C
C     OUTPUT APPROXIMATE SOLUTION X OR A MESSAGE OF FAILURE.
C
      LOGICAL FLAG
      DIMENSION X(3),Z(3),C(3),G(4),A(4)
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG1
      LOGICAL OK
      F1(X1,X2,X3)=3*X1-COS(X2*X3)-.5
      F2(X1,X2,X3)=X1*X1-81*(X2+.1)**2+SIN(X3)+1.06
      F3(X1,X2,X3)=EXP(-X1*X2)+20*X3+(10*PI-3)/3
      F(X1,X2,X3)=F1(X1,X2,X3)**2+F2(X1,X2,X3)**2+F3(X1,X2,X3)**2
C     THE FUNCTION F USED HERE CORRESPONDS TO G USED IN THE ALGORITHM
      PI=3.141593
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Steepest Descent Method'
      WRITE(6,*) 'NOTE THAT THE FUNCTIONS ARE VERY COMPLICATED.'
      WRITE(6,*) 'Have the functions F1,...,Fn and F (G) been'
      WRITE(6,*) 'defined at the beginning of the program and.'
      WRITE(6,*) 'has the gradient been implemented at line 118.'
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AAA
      IF(( AAA .EQ. 'Y' ) .OR. ( AAA .EQ. 'y' )) THEN
         OK = .FALSE.
101      IF (OK) GOTO 11
            WRITE(6,*) 'Input the number n of equations.'
            WRITE(6,*) ' '
            READ(5,*) N
            IF (N.GE.2) THEN
               OK = .TRUE.
            ELSE
               WRITE(6,*) 'N must be greater than 1.'
            ENDIF
         GOTO 101
11       OK = .FALSE.
14       IF (OK) GOTO 15
            WRITE(6,*) 'Input tolerance '
            WRITE(6,*) ' '
            READ(5,*) TOL
            IF (TOL.LE.0.0) THEN
               WRITE(6,*) 'Tolerance must be positive '
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 14
15       OK=.FALSE.
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input the maximum number of iterations.'
            WRITE(6,*) ' '
            READ(5,*) M
            IF ( M .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       DO 16 I=1,N
            WRITE(6,*) 'Input initial approximation X(',I,')'
            WRITE(6,*) ' '
            READ(5,*) X(I)
16       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'can be created '
         OK = .FALSE.
      ENDIF
      IF(.NOT.OK) GOTO 400
      WRITE(6,*) 'Select output destinations: '
      WRITE(6,*) '1. Screen '
      WRITE(6,*) '2. Text file '
      WRITE(6,*) 'Enter 1 or 2 '
      WRITE(6,*) ' '
      READ(5,*) FLAG1
      IF ( FLAG1 .EQ. 2 ) THEN
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
      WRITE(OUP,*) 'STEEPEST DESCENT METHOD FOR NONLINEAR SYSTEMS'
      WRITE(OUP,*) 'Iter. No.   Vector    G(vector)'
C     STEP 1
      K=0
      WRITE(OUP,2) K,X(1),X(2),X(3)
      K=1
C     STEP 2
100   IF (K.GT.M) GOTO 200
C        STEP 3
         C(1)=F1(X(1),X(2),X(3))
         C(2)=F2(X(1),X(2),X(3))
         C(3)=F3(X(1),X(2),X(3))
C        THE VECTOR C HOLDS F(X) FROM THE EXAMPLE AND G1 HOLDS G(X)
         G1=F(X(1),X(2),X(3))
         G(1)=G1
         Q=EXP(-X(1)*X(2))
         P=SIN(X(2)*X(3))
         Z(1)=6*C(1)+4*X(1)*C(2)-2*X(2)*Q*C(3)
         Z(2)=2*X(3)*P*C(1)-324*(X(2)+.1)*C(2)-2*X(1)*Q*C(3)
         Z(3)=2*X(2)*P*C(1)+2*COS(X(3))*C(2)+40*C(3)
         Z0=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
C        STEP 4
         IF(ABS(Z0).LT.1.0E-20) THEN
            WRITE(OUP,8)
            GOTO 400
         ENDIF
C        STEP 5
         DO 10 J=1,N
10       Z(J)=Z(J)/Z0
C        A IS USED FOR ALPHA AND G0 FOR G(3)
         A(1)=0
         X0=1
         G0=F(X(1)-X0*Z(1),X(2)-X0*Z(2),X(3)-X0*Z(3))
C        STEP 6
         FLAG=.TRUE.
500      IF(.NOT.FLAG) GOTO 600
            IF(ABS(G0).LT.ABS(G1)) FLAG=.FALSE.
C           STEPS 7 AND 8
            X0=X0/2
            IF(X0.LT.1.0E-10) THEN
               WRITE(OUP,3)
               GOTO 400
            ENDIF
            G0=F(X(1)-X0*Z(1),X(2)-X0*Z(2),X(3)-X0*Z(3))
         GOTO 500
600      A(3)=X0
         G(3)=G0
C        STEP 9
         X0=X0/2
         G(2)=F(X(1)-X0*Z(1),X(2)-X0*Z(2),X(3)-X0*Z(3))
         A(2)=X0
C        STEP 10
         H1=(G(2)-G(1))/(A(2)-A(1))
         H2=(G(3)-G(2))/(A(3)-A(2))
         H3=(H2-H1)/(A(3)-A(1))
C        STEP 11
         X0=.5*(A(1)+A(2)-H1/H3)
         G0=F(X(1)-X0*Z(1),X(2)-X0*Z(2),X(3)-X0*Z(3))
C        STEP 12
         A0=X0
         DO 20 J=1,N
           IF(ABS(G(J)).LT.ABS(G0)) THEN
               A0=A(J)
               G0=G(J)
           ENDIF
20       CONTINUE
         IF(ABS(A0).LT.1.0E-10) THEN
            WRITE(OUP,9)
            GOTO 400
         ENDIF
C        STEP 13
         DO 30 J=1,N
30       X(J)=X(J)-A0*Z(J)
         WRITE(OUP,2) K,X(1),X(2),X(3),G0
C        STEP 14
         IF(ABS(G0).LT.TOL.OR.ABS(G0-G1).LT.TOL) THEN
            WRITE(OUP,5)
            GOTO 400
         ENDIF
C        STEP 15
         K=K+1
      GOTO 100
C     STEP 16
200   WRITE(OUP,7)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
2     FORMAT(1X,I3,4(1X,E15.8))
8     FORMAT(1X,'ZERO GRADIENT: MAY HAVE MIN')
3     FORMAT(1X,'NO LIKELY IMPROVEMENT: MAY HAVE MIN')
9     FORMAT(1X,'NO CHANGE LIKELY: PROBABLY ROUNDING DIFFICULTIES')
5     FORMAT(1X,'PROCEDURE COMPLETED SUCCESSFULLY')
7     FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED')
      END
