C***********************************************************************
C                                                                      *
C          PIECEWISE LINEAR RAYLEIGH-RITZ METHOD                       *
C                                                                      *
C***********************************************************************
C
C
C
C
C     TO APPROXIMATE THE SOLUTION TO THE BOUNDARY-VALUE PROBLEM
C
C          -D(P(X)Y')/DX + Q(X)Y = F(X), 0 <= X <= 1,
C                  Y(0) = Y(1) = 0
C
C     WITH A PIECEWISE LINEAR FUNCTION:
C
C     INPUT   INTEGER N; MESH POINTS X(0)=0<X(1)<...<X(N)<X(N+1)=1
C
C     OUTPUT  COEFFICIENTS C(1),...,C(N) OF THE BASIS FUNCTIONS
C
C     FUNCTIONS P,Q,F DEFINED AS P,Q,FF IN SUBPROGRAMS
      DIMENSION X(22),H(21),Q(6,21)
      DIMENSION ALPHA(20),BETA(20),B(20),A(20),Z(20),ZETA(20)
      DIMENSION C(20)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is Piecewise Linear Rayleigh-Ritz.'
      WRITE(6,*) 'Have the functions P, QQ and F been created?'
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
14       IF (OK) GOTO 15
            WRITE(6,*) 'Input integer N where X(0)=0, X(N+1)=1'
            WRITE(6,*) ' '
            READ(5,*) N
            IF ( N .LT. 1 ) THEN
              WRITE(6,*) 'N must be greater than 1'
              WRITE(6,*) ' '
            ELSE
              OK = .TRUE.
              N1=N+1
            ENDIF
         GOTO 14
15       OK = .FALSE.
11       IF ( .NOT. OK) THEN
            WRITE(6,*) 'Choice of method to input X(1),...,X(N):'
            WRITE(6,*) '1. Input entry by entry from keyboard '
            WRITE(6,*) '2. Equally spaced nodes to be calculated'
            WRITE(6,*) '3. Input data from a text file '
            WRITE(6,*) 'Choose 1, 2, or 3 please '
            WRITE(6,*) ' '
            READ(5,*)  FLAG
            IF( ( FLAG .GE. 1 ) .AND. ( FLAG .LE. 3 )) OK = .TRUE.
            GOTO 11
         ENDIF
         IF (FLAG .EQ. 1) THEN
            DO 31 I = 2, N1
               J=I-1
               WRITE(6,*) 'Input X(',J,')'
               WRITE(6,*) ' '
               READ(5,*) X(I)
31          CONTINUE
         ENDIF
         IF (FLAG .EQ. 2) THEN
            HC=1.0/(N+1.0)
            DO 21 I=1,N
21          X(I+1) = I*HC
         ENDIF
         IF (FLAG .EQ. 3) THEN
            WRITE(6,*) 'Has a text file been created?'
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
               OK = .TRUE.
               READ(4,*) (X(I),I=2,N1)
               CLOSE(UNIT=4)
            ELSE
               WRITE(6,*) 'The program will end so the input file can '
               WRITE(6,*) 'be created. '
               OK = .FALSE.
            ENDIF
         ENDIF
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'P, QQ and F can be created '
         OK = .FALSE.
      ENDIF
      IF (.NOT.OK) GOTO 400
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
      WRITE(OUP,*) 'PIECEWISE LINEAR RAYLEIGH-RITZ METHOD'
C     SUBSCRIPTS ARE SHIFTED BY ONE TO AVOID ZERO SUBSCRIPTS
      N1 = N+1
      N2 = N+2
      NN = N-1
      X(1)=0.
      X(N2)=1.
      WRITE(OUP,4)
      WRITE(OUP,3) (X(I),I=1,N2)
C     STEP 1
      DO 10 I=1,N1
           H(I)=X(I+1)-X(I)
10    CONTINUE
C     STEP 2
C     PEICEWISE LINEAR BASIS PHI(I) DEFINED IN SUBPROGRAMS AS NEEDED
C     STEP 3
C     COMPUTING THE INTEGRALS FOR THE ENTRIES OF THE MATRIX A
C     SIMPSON'S COMPOSITE METHOD IS USED WITHIN THE SUBPROGRAMS
      DO 20 J=2,N
         Q(1,J-1) = SIMP(1,X(J),X(J+1))/(H(J)*H(J))
         Q(2,J-1) = SIMP(2,X(J-1),X(J))/(H(J-1)*H(J-1))
         Q(3,J-1) = SIMP(3,X(J),X(J+1))/(H(J)*H(J))
         Q(4,J-1) = SIMP(4,X(J-1),X(J))/(H(J-1)*H(J-1))
         Q(5,J-1) = SIMP(5,X(J-1),X(J))/H(J-1)
         Q(6,J-1) = SIMP(6,X(J),X(J+1))/H(J)
20    CONTINUE
      Q(2,N) = SIMP(2,X(N),X(N+1))/(H(N)*H(N))
      Q(3,N) = SIMP(3,X(N+1),X(N+2))/(H(N+1)*H(N+1))
      Q(4,N) = SIMP(4,X(N),X(N+1))/(H(N)*H(N))
      Q(4,N+1) = SIMP(4,X(N+1),X(N+2))/(H(N+1)*H(N+1))
      Q(5,N) = SIMP(5,X(N),X(N+1))/H(N)
      Q(6,N) = SIMP(6,X(N+1),X(N+2))/H(N+1)
C     STEP 4
      DO 30 J = 1,NN
         ALPHA(J) = Q(4,J)+Q(4,J+1)+Q(2,J)+Q(3,J)
         BETA(J) = Q(1,J)-Q(4,J+1)
         B(J) = Q(5,J)+Q(6,J)
30    CONTINUE
C     STEP 5
      ALPHA(N) = Q(4,N)+Q(4,N+1)+Q(2,N)+Q(3,N)
      B(N) = Q(5,N)+Q(6,N)
C     STEP 6
      A(1) = ALPHA(1)
      ZETA(1) = BETA(1)/ALPHA(1)
C     STEP 7
      DO 40 J = 2,NN
         A(J) = ALPHA(J)-BETA(J-1)*ZETA(J-1)
         ZETA(J) = BETA(J)/A(J)
40    CONTINUE
C     STEP 8
      A(N) = ALPHA(N)-BETA(N-1)*ZETA(N-1)
C     STEP 9
      Z(1) = B(1)/A(1)
C     STEP 10
      DO 50 J=2,N
         Z(J) = (B(J)-BETA(J-1)*Z(J-1))/A(J)
50    CONTINUE
C     STEP 11
      C(N) = Z(N)
C     STEP 12
      DO 60 J=1,NN
         J1 = N-J
         C(J1) = Z(J1)-ZETA(J1)*C(J1+1)
60    CONTINUE
      WRITE(OUP,5)
      WRITE(OUP,3) (C(I),I=1,N)
C     STEP 13
C     PROCESS IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT (8(1X,E15.8))
5     FORMAT (1X,'OUTPUT IS C(1), C(2), ... , C(N)')
4     FORMAT (1X,'OUTPUT IS X(0), X(1), ... , X(N), X(N+1)')
      END
      FUNCTION P(X)
      P = 1.0
      RETURN
      END
      FUNCTION QQ(X)
      QQ = (3.141593)**2
      RETURN
      END
      FUNCTION FF(X)
      PI = 3.141593
      FF = 2.0*PI*PI*SIN(PI*X)
      RETURN
      END
      FUNCTION SIMP(NF,A,B)
      DIMENSION Z(5)
      H = (B-A)/4.0
      DO 10 I=1,5
         Y = A+I*H
         IF(NF.EQ.1) Z(I) = (4.0-I)*I*H*H*QQ(Y)
         IF(NF.EQ.2) Z(I) = I*I*H*H*QQ(Y)
         IF(NF.EQ.3) Z(I) = (H*(4.0-I))**2*QQ(Y)
         IF(NF.EQ.4) Z(I) = P(Y)
         IF(NF.EQ.5) Z(I) = I*H*FF(Y)
         IF(NF.EQ.6) Z(I) = (4.0-I)*H*FF(Y)
10    CONTINUE
      SIMP = (Z(1)+Z(5)+2*Z(3)+4*(Z(2)+Z(4)))*H/3
      RETURN
      END
