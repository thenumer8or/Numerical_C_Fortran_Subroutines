C***********************************************************************
C                                                                      *
C                        SOR METHOD                                    *
C                                                                      *
C***********************************************************************
C
C
C
C     TO SOLVE AX = B GIVEN THE PARAMETER W AND AN ITITIAL APPROXIMATION
C     X(O):
C
C     INPUT:   THE NUMBER OF EQUATIONS AND UNKNOWNS N; THE ENTRIES
C              A(I,J), 1<=I, J<=N, OF THE MATRIX A; THE ENTRIES B(I),
C              1<=I<=N, OF THE INHOMOGENEOUS TERM B; THE ENTRIES XO(I),
C              1<=I<=N, OF X(0); THE PARAMETER W; TOLERANCE TOL;
C              MAXIMUM NUMBER OF ITERATIONS N.
C
C     OUTPUT:  THE APPROXIMATE SOLUTION X(1),...,X(N) OR A MESSAGE
C              THAT THE NUMBER OF ITERATIONS WAS EXCEEDED.
C
C     INITIALIZATION
      DIMENSION A(3,4),X1(3)
C     USE NN FOR CAPITAL N
C     USE W FOR OMEGA
C     B(I) = A(I,N+1) FOR 1<=I<=N
C     USE X1 FOR XO
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the SOR Method for Linear Systems.'
      WRITE(6,*) 'The array will be input from a text file in the'
      WRITE(6,*) ' order: A(1,1), A(1,2), ..., A(1,N+1), A(2,1),'
      WRITE(6,*) ' A(2,2), ..., A(2,N+1)..., A(N,1), A(N,2),'
      WRITE(6,*) ' ..., A(N,N+1) '
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
         WRITE(6,*) 'Input the number of equations - an integer '
         WRITE(6,*)
         READ(5,*) N
         IF (N .GT. 0) THEN
            M = N+1
            READ(INP,*) ((A(I,J), J=1,M),I=1,N)
            READ(INP,*) (X1(I),I=1,N)
          OK = .TRUE.
            CLOSE(UNIT=INP)
         ELSE
            WRITE(6,*) 'The number must be a positive integer'
         ENDIF
         GOTO 19
111      OK = .FALSE.
12       IF (OK) GOTO 13
         WRITE(6,*) 'Input the tolerance.'
         WRITE(6,*) ' '
         READ(5,*) TOL
         IF (TOL .GT. 0.0) THEN
            OK = .TRUE.
         ELSE
            WRITE(6,*) 'Tolerance must be positive.'
         ENDIF
         GOTO 12
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
15       WRITE(6,*) 'Input parameter w (omega).'
         WRITE(6,*) ' '
         READ(5,*) W
         GOTO 16
      ELSE
         WRITE(6,*) 'The program will end so the input file can '
         WRITE(6,*) 'be created. '
         OK = .FALSE.
      ENDIF
16    IF(.NOT. OK) GOTO 400
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
      WRITE(OUP,*) 'SOR ITERATIVE METHOD FOR LINEAR SYSTEMS'
      WRITE(OUP,3)
      WRITE(OUP,4) ((A(I,J),J=1,M),I=1,N)
      WRITE(OUP,5)
      WRITE(OUP,4) ( X1(I) ,I=1,N)
      WRITE(OUP,11) W
C     STEP 1
      K=1
C     STEP 2
100   IF (K.GT.NN) GOTO 200
           ERR=0.0
C          ERR WILL BE USED TO TEST ACCURACY, IT MEASURES THE L2 NORM
C          STEP 3
C          THE DO-LOOP COMPUTES THE SUMMATION
           DO 10 I=1,N
                S=0.0
                DO 20 J=1,N
20                   S=S-A(I,J)*X1(J)
                S=W*(S+A(I,N+1))/A(I,I)
                ERR=ERR+S*S
C               X IS NOT USED, SINCE ONLY ONE VECTOR IS NEEDED
10              X1(I)=X1(I)+S
           ERR=SQRT(ERR)
           WRITE(OUP,6) K,ERR,(X1(I),I=1,N)
C          STEP 4
           IF(ERR.LE.TOL) THEN
                WRITE(OUP,7) K
                GOTO 400
           END IF
C          STEP 5
           K=K+1
C          STEP 6 IS NOT NEEDED
      GOTO 100
C     STEP 7
C     PROCEDURE COMPLETED SUCCESSFULLY
200   CONTINUE
      WRITE(OUP,8) NN
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'ORIGINAL SYSTEM:')
4     FORMAT((1X,4(3X,E15.8)))
5     FORMAT(1X,'INITIAL APPROXIMATION:')
6     FORMAT(1X,'ITERATION NUMBER ',I3,' GIVES ERROR ',E15.8,' FOR APPRO
     *X. ',3(2X,E15.8))
7     FORMAT(1X,'CONVERGENCE ON INTERATION NUMBER ',I4)
8     FORMAT(1X,'FAILURE AFTER INTERATION NUMBER ',I4)
11    FORMAT(1X,'OMEGA IS ',E15.8)
      END
