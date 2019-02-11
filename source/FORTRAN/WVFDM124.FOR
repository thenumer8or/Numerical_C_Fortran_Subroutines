C***********************************************************************
C                                                                      *
C             WAVE EQUATION FINITE-DIFFERENCE METHOD                   *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION TO THE WAVE EQUTION:
C     SUBJECT TO THE BOUNDARY CONDITIONS
C              U(0,T) = U(L,T) = 0, 0 < T < MAX T
C     AND THE INITIAL CONDITIONS
C          U(X,0) = F(X) AND DU(X,0)/DT = G(X), 0 <= X <= L:
C
C     INPUT:   ENDPOINT L; MAXIMUM TIME T; CONSTANT ALPHA; INTEGERS M,N.
C
C     OUTPUT:  APPROXIMATIONS W(I,J) TO U(X(I),T(J)) FOR EACH I=0,...M
C              AND J=0,..,N.
C
C     INITIALIZATION
C     DIMENSION W(M+1,N+1)
      DIMENSION W(11,21)
C     FT IS CAPITAL T AND FX IS L
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      F(XZ)=SIN(PI*XZ)
      G(XZ)=0
      PI=3.141593
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Finite-Difference Method'
      WRITE(6,*) 'for the Wave Equation.'
      WRITE(6,*) 'Have the functions F and G been created?'
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
         WRITE(6,*) 'The lefthand endpoint on the X-axis is 0.'
10       IF (OK) GOTO 11
            WRITE(6,*) 'Input righthand endpoint on the X-axis.'
            WRITE(6,*) ' '
            READ(5,*) FX
            IF (FX.lE.0.0) THEN
               WRITE(6,*) 'Must be a positive number.'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 10
11       OK = .FALSE.
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input the maximum value of the time'
            WRITE(6,*) 'variable T.'
            WRITE(6,*) ' '
            READ(5,*) FT
            IF ( FT .LE. 0.0 ) THEN
              WRITE(6,*) 'Must be positive number.'
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       WRITE(6,*) 'Input the constant alpha.'
         WRITE(6,*) ' '
         READ(5,*) ALPHA
15       OK = .FALSE.
14       IF (OK) GOTO 16
            WRITE(6,*) 'Input integer m = number of intervals'
            WRITE(6,*) 'on X-axis and N = number of time intervals.'
            WRITE(6,*) 'Note: m should exceed 3. Separate by blank.'
            WRITE(6,*) ' '
            READ(5,*) M,N
            IF ((N.LE.0).OR.(N.LE.0)) THEN
              WRITE(6,*) 'Must be positive integers '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 14
16       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'F and G can be created.'
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
      WRITE(OUP,*) 'FINITE DIFFERENCE METHOD FOR WAVE EQUATION'
      M1=M+1
      N1=N+1
      MM=M-1
      NN=N-1
C     STEP 1
C     TK IS USED FOR K, V FOR LAMBDA
      H=FX/M
      TK=FT/N
      V=TK*ALPHA/H
C     STEP 2
C     THE SUBSCRIPTS ARE SHIFTED TO AVOID ZERO SUBSCRIPTS
      DO 20 J=2,N1
           W(1,J)=0
20    W(M1,J)=0
C     STEP 3
      W(1,1)=F(0.0)
      W(M1,1)=F(FX)
C     STEP 4
      DO 30 I=2,M
           W(I,1)=F((I-1)*H)
30    W(I,2)=(1-V*V)*F((I-1)*H)+V*V*(F(I*H)+F((I-2)*H))/2+TK*G((I-1)*H)
C     STEP 5
      DO 40 J=2,N
           DO 40 I=2,M
40    W(I,J+1)=2*(1-V*V)*W(I,J)+V*V*(W(I+1,J)+W(I-1,J))-W(I,J-1)
C     STEP 6
C     OUTPUT ONLY FOR TIME VALUE T=FT
      WRITE(OUP,1) FT
      WRITE(OUP,3)
      DO 50 I=1,M1
           X=(I-1)*H
50    WRITE(OUP,2) I,X,W(I,N1)
C     STEP 7
C     PROCEDURE IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'FINAL TIME VALUE IS',1X,E15.8)
2     FORMAT(1X,I3,2(1X,E15.8))
3     FORMAT(3X,'I',12X,'T(I)',10X,'W(I,N)')
      END
