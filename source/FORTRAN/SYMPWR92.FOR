C***********************************************************************
C                                                                      *
C           SYMMETRIC POWER METHOD                                     *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE DOMINANT EIGENVALUE AND AN ASSOCIATED
C     EIGENVECTOR OF THE N BY N SYMMETRIC MATRIX A GIVEN A NONZERO:
C     VECTOR X:
C
C     INPUT:   DIMENSION N; MATRIX A; VECTOR X; TOLERANCE TOL; MAXIMUM
C              NUMBER OF ITERATIONS N.
C
C     OUTPUT:  APPROXIMATE EIGENVALUE MU; APPROXIMATE EIGENVECTOR X OR
C              A MESSAGE THAT THE MAXIMUM NUMBER OF ITERATIONS WAS
C              EXCEEDED.
C
C     INITIALIZATION
      DIMENSION A(10,10),X(10),Y(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Symmetric Power Method.'
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
            READ(INP,*) (Y(I),I=1,N)
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
      WRITE(OUP,*) 'SYMMETRIC POWER METHOD'
      DO 30 I=1,N
30    X(I) = 0.0
C     USE NN IN PLACE OF N FOR THE MAXIMUM NUMBER OF ITERATIONS
      WRITE(OUP,3)
      WRITE(OUP,4) ((A(I,J),J=1,N),I=1,N)
      WRITE(OUP,5)
      WRITE(OUP,4) (Y(I),I=1,N)
C     STEP 1
      K = 1
C     NORM COMPUTES THE NORM OF THE VECTOR X AND RETURNS X DIVIDED BY
C     ITS NORM IN THE VARIABLE Y
      CALL NORM(N,Y,X,ERR)
C     STEP 2
100   IF (K.GT.NN) GOTO 200
C     STEPS 3 AND 4
           YMU=0
           DO 10 I=1,N
                Y(I) = 0.0
                DO 20 J=1,N
20                   Y(I) = Y(I)+A(I,J)*X(J)
10              YMU = YMU+X(I)*Y(I)
C          STEP 5 THIS STEP IS ACCOMPLISHED IN SUBROUTINE NORM
C          STEP 6
           CALL NORM(N,Y,X,ERR)
           WRITE(OUP,7) K,YMU
           WRITE(OUP,8) (Y(I),I=1,N)
           WRITE(OUP,9) (X(I),I=1,N)
C          STEP 7
           IF(ERR.LT.TOL) THEN
                WRITE(OUP,6) YMU
C               PROCEDURE COMPLETED SUCCESSFULLY
                GOTO 400
           END IF
C          STEP 8
           K = K+1
      GOTO 100
C     STEP 9
200   WRITE(OUP,11) NN
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'MATRIX A')
4     FORMAT(1X,3(5X,E15.8))
5     FORMAT(1X,'INITIAL VECTOR')
6     FORMAT(1X,'EIGENVALUE IS ',E15.8)
7     FORMAT(1X,'ITERATION NUM. ',I3,' GIVES APPROX.= ',E15.8)
8     FORMAT(1X,'Y VECTOR IS',4(5X,E15.8))
9     FORMAT(1X,'UNIT X VECTOR IS',4(5X,E15.8))
11    FORMAT(1X,'FAILURE AFTER ITERATION NUM. ',I3)
      END
C
           SUBROUTINE NORM(N,Y,X,ERR)
           DIMENSION X(N),Y(N)
           XL = 0.0
           DO 10 I=1,N
10         XL = XL+Y(I)*Y(I)
           XL = SQRT(XL)
           ERR=0.0
           IF(XL.GE.1.0E-20) THEN
                DO 20 I=1,N
                   T=Y(I)/XL
                   IF(ABS(X(I)-T).GT.ERR) ERR=ABS(X(I)-T)
20              X(I) = T
                RETURN
           END IF
           WRITE(6,1)
1          FORMAT(1X,'ERROR-Y IS ZERO, SO A HAS ZERO EIGENVALUE')
           STOP
           END
