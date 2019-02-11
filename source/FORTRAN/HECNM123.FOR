C***********************************************************************
C                                                                      *
C                    CRANK-NICOLSON METHOD                             *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION TO A PARABOLIC PARTIAL-DIFFERENTIAL
C     EQUATION SUBJECT TO THE BOUNDARY CONDITIONS
C                U(0,T) = U(L,T) = 0, 0 < T < MAX T
C     AND THE INITIAL CONDITIONS
C                 U(X,0) = F(X), 0 <= X <= L:
C
C     INPUT:   ENDPOINT L; MAXIMUM TIME T; CONSTANT ALPHA; INTEGERS
C              M,N:
C
C     OUTPUT:  APPROXIMATIONS W(I,J) TO U(X(I),T(J)) FOR EACH
C              I=1,...,M-1 AND J=1,...,N.
C
C     INITIALIZATION
C     DIMENSION V(M),XL(M-1),XU(M-1),Z(M-1)
      DIMENSION V( 10),XL( 9),XU( 9),Z( 9)
C     V IS USED FOR W, FT FOR CAPITAL T, FX FOR L
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      F(XZ)=SIN(PI*XZ)
      PI=3.1415927
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Crank-Nicolson Method'
      WRITE(6,*) 'for the Heat Equation.'
      WRITE(6,*) 'Has the function F been created in the program? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
         WRITE(6,*) 'The lefthand endpoint on the X-axis is 0.'
19       IF (OK) GOTO 11
            WRITE(6,*) 'Input righthand endpoint on the X-axis.'
            WRITE(6,*) ' '
            READ(5,*) FX
            IF (FX.lE.0.0) THEN
               WRITE(6,*) 'Must be a positive number.'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 19
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
         WRITE(6,*) 'The program will end so that the function F '
         WRITE(6,*) 'can be created '
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
      WRITE(OUP,*) 'CRANK-NICOLSON METHOD FOR HEAT EQUATION'
      MM=M-1
      MMM=MM-1
      NN=N-1
C     STEP 1
      H=FX/M
C     TK IS USED FOR K
      TK=FT/N
C     VV IS USED FOR LAMBDA
      VV=ALPHA*ALPHA*TK/(H*H)
C     SET V(M)=0
      V(M)=0
C     STEP2
C     INITIAL VALUES
      DO 10 I=1,MM
10         V(I)=F(I*H)
C     STEP 3
C     STEPS 3 THROUGH 11 SOLVE A TRIDIAGONAL LINEAR SYSTEM
C     USING ALGORITHM 6.7
C     USE XL FOR L, XU FOR U
      XL(1)=1+VV
      XU(1)=-VV/(2*XL(1))
C     STEP 4
      DO 20 I=2,MMM
           XL(I)=1+VV+VV*XU(I-1)/2
20    XU(I)=-VV/(2*XL(I))
C     STEP 5
      XL(MM)=1+VV+VV*XU(MMM)/2
C     STEP 6
      DO 30 J=1,N
C          STEP 7
C          CURRENT T(J)
           T=J*TK
           Z(1)=((1-VV)*V(1)+VV*V(2)/2)/XL(1)
C          STEP 8
           DO 40 I=2,MM
40         Z(I)=((1-VV)*V(I)+VV/2*(V(I+1)+V(I-1)+Z(I-1)))/XL(I)
C          STEP 9
           V(MM)=Z(MM)
C          STEP 10
           DO 50 II=1,MMM
                I=MMM-II+1
50         V(I)=Z(I) -XU(I)*V(I+1)
30    CONTINUE
C     STEP 11--OUTPUT WILL BE ONLY FOR T=FT
      WRITE(OUP,1) FT
      WRITE(OUP,3)
      DO 60 I=1,MM
           X=I*H
60    WRITE(OUP,2) I,X,V(I)
C     STEP 12
C     PROCEDURE IS COMPLETE
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.OUP) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'OUTPUT FOR TIME = ' ,1X,E15.8)
2     FORMAT(1X,I3,2(1X,E15.8))
3     FORMAT(3X,'I',12X,'T(I)',12X,'V(I)')
      END
