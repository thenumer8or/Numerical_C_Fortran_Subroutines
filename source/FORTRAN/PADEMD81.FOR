C***********************************************************************
C                                                                      *
C        PADE RATIONAL APPROXIMATION METHOD                            *
C                                                                      *
C***********************************************************************
C
C
C
C     To obtain the rational approximation
C
C      r(x) = p(x) / q(x)
C           = (p0 + p1*x + ... + Pn*x^n) / (q0 + q1*x + ... + qm*x^m)
C
C     for a given function f(x):
C
C     INPUT  nonnegative integers m and n.
C
C     OUTPUT  coefficients q0, q1, ... , qm, p0, p1, ... , pn.
C
      DIMENSION A(10,11),AA(11),NROW(10),P(7),Q(7)
      INTEGER PP,BN
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      ZERO = 1.0E-20
      OK = .FALSE.
      WRITE(6,*) 'This is Pade Approximation.'
11    IF (OK) GOTO 12
      WRITE(6,*) 'Input m and n'
      WRITE(6,*) ' '
      READ(5,*) LM,LN
      BN = LM+LN
      IF( (LM.GE.0) .AND. (LN.GE.0)) THEN
         OK = .TRUE.
      ELSE
         WRITE(6,*) 'm and n must be nonnegative'
      ENDIF
      GOTO 11
12    IF ( (LM.EQ.0) .AND. (LN.EQ.0)) THEN
         OK = .FALSE.
         WRITE(6,*) 'Not both m and n can be zero.'
         GOTO 11
      ENDIF
      OK = .FALSE.
21    IF ( .NOT. OK) THEN
         WRITE(6,*) 'The McLaurin coefficients a(0),...,a(N)'
         WRITE(6,*) 'are to be input.'
         WRITE(6,*) 'Choice of input method: '
         WRITE(6,*) '1. Input entry by entry from keyboard '
         WRITE(6,*) '2. Input data from a text file '
         WRITE(6,*) 'Choose 1 or 2 please '
         WRITE(6,*) ' '
         READ(5,*)  FLAG
         IF( ( FLAG .GE. 1 ) .AND. ( FLAG .LE. 2 )) OK = .TRUE.
         GOTO 21
      ENDIF
      IF (FLAG .EQ. 1) THEN
         DO 31 I = 0, BN
            WRITE(6,*) 'Input a(',I,')'
            WRITE(6,*) ' '
            READ(5,*) AA(I+1)
31       CONTINUE
      ENDIF
      IF (FLAG .EQ. 2) THEN
         WRITE(6,*) 'Has a text file been created.'
         WRITE(6,*) 'Enter Y or N - letter within quotes '
         WRITE(6,*) ' '
         READ(5,*)  AAA
         IF (( AAA .EQ. 'Y' ) .OR.( AAA .EQ. 'y' )) THEN
            WRITE(6,*) 'Input the file name in the form - '
            WRITE(6,*) 'drive:name.ext  contained in quotes'
            WRITE(6,*) 'as example:   ''A:DATA.DTA'' '
            WRITE(6,*) ' '
            READ(5,*)  NAME
            INP = 4
            OPEN(UNIT=INP,FILE=NAME,ACCESS='SEQUENTIAL')
            READ(4,*) (AA(I+1),I=0,BN)
            CLOSE(UNIT=4)
         ELSE
            WRITE(6,*) 'The program will end so the input file can '
            WRITE(6,*) 'be created. '
            OK = .FALSE.
         ENDIF
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
      WRITE(OUP,*) 'PADE APPROXIMATION'
C     STEP 1
      N = BN
      M = N + 1
C     STEP 2  -  performed in input
      DO 20 I = 1, N
20       NROW(I) = I
C     initialize row pointer for linear system
      NN = N - 1
C     STEP 3
      Q(1) = 1.0
      P(1) = AA(1)
C     STEP 4
C     Set up a linear system, but use A[i,j] instead of B[i,j].
      DO 30 I = 1, N
C        STEP 5
         LLL = I - 1
         IF (LLL .GE. 1) THEN
         DO 40 J = 1, LLL
            IF (J .LE. LN) A(I,J) = 0.0
40       CONTINUE
      ENDIF
C        STEP 6
         IF (I .LE. LN) A(I,I) = 1.0
C        STEP 7
         LLL = I + 1
         IF (LLL .LE. LN) THEN
         DO 50 J = LLL, LN
50          A(I,J) = 0.0
         ENDIF
C        STEP 8
         IF (I .GE. 1 ) THEN
         DO 60 J = 1, I
60          IF (J .LE. LM) A(I,LN+J) = -AA(I-J+1)
         ENDIF
C        STEP 9
         LLL = LN+I+1
         IF (LLL .LE. N) THEN
         DO 70 J = LLL, N
70          A(I,J) = 0.0
         ENDIF
C        STEP 10
         A(I,N+1) = AA(I+1)
30    CONTINUE
      ICHG = 0
C     Solve the linear system using partial pivoting.
      I = LN+1
C     STEP 11
80    IF ((.NOT. OK ) .OR. ( I .GT. NN )) GOTO 90
C        STEP 12
         IMAX = NROW(I)
         AMAX = ABS( A(IMAX,I) )
         IMAX = I
         JJ = I + 1
         DO 100 IP = JJ, N
            JP = NROW(IP)
            IF ( ABS( A(JP,I) ) .GT. AMAX ) THEN
               AMAX = ABS( A(JP,I) )
               IMAX = IP
            ENDIF
100      CONTINUE
C        STEP 13
         IF ( AMAX .LE. ZERO ) THEN
            OK = .FALSE.
         ELSE
C           STEP 14
C           simulate row interchange
            IF ( NROW(I) .NE. NROW(IMAX) ) THEN
               ICHG = ICHG + 1
               NCOPY = NROW(I)
               NROW(I) = NROW(IMAX)
               NROW(IMAX) = NCOPY
            ENDIF
            I1 = NROW(I)
C           STEP 15
C           Perform elimination.
            DO 110 J = JJ, M
               J1 = NROW(J)
C              STEP 16
               XM = A(J1,I) / A(I1,I)
C              STEP 17
               DO 120 K = JJ, M
120               A(J1,K) = A(J1,K) - XM * A(I1,K)
C              STEP 18
               A(J1,I) = XM
110         CONTINUE
         ENDIF
         I = I + 1
      GOTO 80
90    IF ( OK ) THEN
C        STEP 19
         N1 = NROW(N)
         IF ( ABS( A(N1,N) ) .LE. ZERO ) THEN
            OK = .FALSE.
C           system has no unique solution
         ELSE
C           STEP 20
C           Start backward substitution.
            IF (LM .GT. 0) THEN
               Q(LM+1) = A(N1,M) / A(N1,N)
               A(N1,M) = Q(LM+1)
            ENDIF
            PP = 1
C           STEP 21
            LLL = LN+1
            IF (LLL .LE. NN) THEN
            DO 130 K = LLL, NN
               I = NN - K + LN+1
               JJ = I + 1
               N2 = NROW(I)
               SUM = A(N2,N+1)
               IF (JJ .LE. N) THEN
               DO 140 KK = JJ, N
                  LL = NROW(KK)
                  SUM = SUM - A(N2,KK) * A(LL,M)
140            CONTINUE
               ENDIF
               A(N2,M) = SUM / A(N2,I)
               Q(LM-PP+1) = A(N2,M)
               PP = PP + 1
130         CONTINUE
            ENDIF
C           STEP 22
            DO 150 K = 1, LN
               I = LN - K + 1
               N2 = NROW(I)
               SUM = A(N2,N+1)
               LLL = LN+1
               IF (LLL .LE. N) THEN
               DO 160 KK = LLL, N
                  LL = NROW(KK)
                  SUM = SUM - A(N2,KK) * A(LL,M)
160            CONTINUE
               ENDIF
               A(N2,M) = SUM
               P(LN-K+2) = A(N2,M)
150         CONTINUE
C           STEP 23
C           procedure completed successfully
         ENDIF
      ENDIF
      IF ( .NOT. OK ) THEN
         WRITE(OUP,*) 'System has no unique solution '
      ELSE
         WRITE(OUP,*)  'Coefficients Q[0], ..., Q[M]'
         DO 170 I = 0, LM
170         WRITE(OUP,1) Q(I+1)
         WRITE(OUP,*)  'Coefficients P[0], ..., P[N]'
         DO 180 I = 0, LN
180         WRITE(OUP,1) P(I+1)
      ENDIF
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,E15.8)
      END
