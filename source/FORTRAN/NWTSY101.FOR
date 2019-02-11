C***********************************************************************
C                                                                      *
C             NEWTON'S METHOD FOR SYSTEMS                              *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION OF THE NONLINEAR SYSTEM F(X)=0 GIVEN
C     AN INITIAL APPROXIMATION X:
C
C     INPUT:   NUMBER N OF EQUATIONS AND UNKNOWNS; INITIAL APPROXIMATION
C              X=(X(1),...,X(N)); TOLERANCE TOL; MAXIMUM NUMBER OF
C              ITERATIONS N.
C
C     OUTPUT:  APPROXIMATE SOLUTION X=(X(1),...,X(N)) OR A MESSAGE
C              THAT THE NUMBER OF ITERATIONS WAS EXCEEDED.
C
C     INITIALIZATION
      DIMENSION AA(3,4),Y(3),X(3)
C     USE AA FOR J;  NN FOR MAXIMUM NUMBER OF ITERATIONS N
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Newton Method for Nonlinear Systems.'
      WRITE(6,*) 'Has the function F been defined and has the'
      WRITE(6,*) 'Jacobian of partial derivatives been defined'
      WRITE(6,*) 'at lines 104 and 114 respectively.'
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
            READ(5,*) NN
            IF ( NN .LE. 0 ) THEN
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
         WRITE(6,*) 'F(X) and J(X) can be created '
         OK = .FALSE.
      ENDIF
      IF(.NOT.OK) GOTO 400
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
      WRITE(OUP,*) 'NEWTON METHOD FOR NONLINEAR SYSTEMS'
      PI = 4*ATAN(1.)
C     STEP 1
      K=1
C     STEP 2
100   IF (K.GT.NN) GOTO 200
C          STEP3
C          COMPUTE J(X)
           AA(1,1)=3.0
           AA(1,2)=X(3)*SIN(X(2)*X(3))
           AA(1,3)=X(2)*SIN(X(2)*X(3))
           AA(2,1)=2*X(1)
           AA(2,2)=-162*(X(2)+.1)
           AA(2,3)=COS(X(3))
           AA(3,1)=-X(2)*EXP(-X(1)*X(2))
           AA(3,2)=-X(1)*EXP(-X(1)*X(2))
           AA(3,3)=20
C          COMPUTE -F(X)
           AA(1,4)=-(3*X(1)-COS(X(2)*X(3))-.5)
           AA(2,4)=-(X(1)**2-81*(X(2)+.1)**2+SIN(X(3))+1.06)
           AA(3,4)=-(EXP(-X(1)*X(2))+20*X(3)+(10*PI-3)/3)
C          STEP 4
C          SOLVES THE N X N LINEAR SYSTEM
           CALL LIN(N,N+1,AA,Y)
C          STEPS 5 AND 6
C          R = INFINITY NORM OF Y
           R=0
           DO 10 I=1,N
                IF(ABS(Y(I)).GT.R) R=ABS(Y(I))
10              X(I)=X(I)+Y(I)
           WRITE(OUP,3) K,(X(I),I=1,N),R
C          STEP 6
           IF(R.LT.TOL) THEN
C               PROCESS IS COMPLETE
                WRITE(OUP,4)
                GOTO 400
           END IF
C          STEP 7
           K=K+1
      GOTO 100
C     STEP 8
C     DIVERGENCE
200   WRITE(OUP,5) NN
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'ITER.=',1X,I2,1X,'APPROX. SOL. IS',1X,3(1X,E15.8),/,1X,
     *'APPROX. ERROR IS',1X,E15.8)
4     FORMAT(1X,'SUCCESS WITHIN TOLERANCE 1.0E-05')
5     FORMAT(1X,'DIVERGENCE - STOPPED AFTER ITER.',1X,I2)
      END
C
           SUBROUTINE LIN(N,M,A,X)
           DIMENSION A(N,M),X(N)
           K=N-1
           DO 10 I=1,K
                Y=ABS(A(I,I))
                IR=I
                IA=I+1
                DO 20 J=IA,N
                     IF(ABS(A(J,I)).GT.Y) THEN
                          IR=J
                          Y=ABS(A(J,I))
                     END IF
20              CONTINUE
                IF(Y.EQ.0.0) THEN
                     WRITE(6,1)
                     STOP
                END IF
                IF(IR.NE.I) THEN
                     DO 30 J=I,M
                          C=A(I,J)
                          A(I,J)=A(IR,J)
30                   A(IR,J)=C
                END IF
                DO 40 J=IA,N
                     C=A(J,I)/A(I,I)
                     DO 40 L=I,M
                     IF(ABS(C).LT.1.0E-25) C=0
40              A(J,L)=A(J,L)-C*A(I,L)
10         CONTINUE
           IF(ABS(A(N,N)).LT.1.0E-20) THEN
                WRITE(6,1)
                STOP
           END IF
      X(N)=A(N,M)/A(N,N)
      DO 50 I=1,K
           J=N-I
           JA=J+1
           C=A(J,M)
           DO 60 L=JA,N
60        C=C-A(J,L)*X(L)
50    X(J) = C/A(J,J)
      RETURN
1     FORMAT(1X,'LINEAR SYSTEM HAS NO SOLUTION')
      END
