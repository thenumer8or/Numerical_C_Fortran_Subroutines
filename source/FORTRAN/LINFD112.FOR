C***********************************************************************
C                                                                      *
C               LINEAR FINITE-DIFFERENCE METHOD                        *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION OF THE BOUNDARY-VALUE PROBLEM
C
C           Y'' = P(X)Y' + Q(X)Y + R(X), A <= X <= B
C                 Y(A) = ALPHA, Y(B) = BETA:
C
C     INPUT:   ENDPOINTS A,B; BOUNDARY CONDITIONS ALPHA, BETA;
C              INTEGER N.
C
C     OUTPUT:  APPROXIMATIONS W(I) TO Y(X(I)) FOR EACH I=0,1,...N+1.
C
C     INITIALIZATION
      DIMENSION A(9),B(9),C(9),D(9),XL(9),XU(9),Z(9),W(9)
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      P(X)=-2/X
      Q(X)=2/X**2
      R(X)=SIN(ALOG(X))/X**2
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Linear Finite Difference Method.'
      WRITE(6,*) 'Have the functions P, Q and R been created? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AAA
      IF(( AAA .EQ. 'Y' ) .OR. ( AAA .EQ. 'y' )) THEN
         OK = .FALSE.
110      IF (OK) GOTO 11
            WRITE(6,*) 'Input left and right endpoints separated by'
            WRITE(6,*) 'blank - include decimal point'
            WRITE(6,*) ' '
            READ(5,*) AA, BB
            IF (AA.GE.BB) THEN
               WRITE(6,*) 'Left endpoint must be less'
               WRITE(6,*) 'than right endpoint'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 110
11       OK = .FALSE.
         WRITE(6,*) 'Input Y(',AA,')'
         WRITE(6,*) ' '
         READ(5,*) ALPHA
         WRITE(6,*) 'Input Y(',BB,')'
         WRITE(6,*) ' '
         READ(5,*) BETA
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input a positive integer for the number'
            WRITE(6,*) 'of subinvervals. note: h = (b-a)/(n+1) '
            WRITE(6,*) ' '
            READ(5,*) N
            IF ( N .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'P, Q and R can be created '
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
      WRITE(OUP,*) 'LINEAR FINITE DIFFERENCE METHOD'
C     STEP 1
      H=(BB-AA)/(N+1)
      X=AA+H
      A(1)=2+H*H*Q(X)
      B(1)=-1+H*P(X)/2
      D(1)=-H*H*R(X)+(1+H*P(X)/2)*ALPHA
      M=N-1
C     STEP 2
      DO 10 I=2,M
           X=AA+I*H
           A(I)=2+H*H*Q(X)
           B(I)=-1+H*P(X)/2
           C(I)=-1-H*P(X)/2
10    D(I)=-H*H*R(X)
C     STEP 3
      X=BB-H
      A(N)=2+H*H*Q(X)
      C(N)=-1-H*P(X)/2
      D(N)=-H*H*R(X)+(1-H*P(X)/2)*BETA
C     STEP 4
C     STEPS 4 THROUGH 10 SOLVE A TRIADIAGONAL LINEAR SYSTEM USING
C     ALGORITHM 6.7
C     USE XL, XU FOR L, U RESP.
      XL(1)=A(1)
      XU(1)=B(1)/A(1)
C     STEP 5
      DO 20 I=2,M
           XL(I)=A(I)-C(I)*XU(I-1)
20    XU(I)=B(I)/XL(I)
C     STEP 6
      XL(N)=A(N)-C(N)*XU(N-1)
C     STEP 7
      Z(1)=D(1)/XL(1)
C     STEP 8
      DO 30 I=2,N
30    Z(I)=(D(I)-C(I)*Z(I-1))/XL(I)
C     STEP 9
      W(N)=Z(N)
C     STEP 10
      DO 40 J=1,M
           I=N-J
40    W(I)=Z(I)-XU(I)*W(I+1)
C     STEP 11
      WRITE(OUP,1)
      WRITE(OUP,2) AA,ALPHA
      DO 50 I=1,N
           X=AA+I*H
50    WRITE(OUP,2) X,W(I)
      WRITE(OUP,2) BB,BETA
C     STEP 12
C     PROCEDURE IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'ORDER OF OUTPUT - X(I),W(I)')
2     FORMAT(1X,2(E15.8,3X))
      END
