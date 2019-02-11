C***********************************************************************
C                                                                      *
C   CROUT REDUCTION FOR TRIDIAGONAL LINEAR SYSTEMS                     *
C                                                                      *
C***********************************************************************
C
C
C
C     TO SOLVE AN N BY N LINEAR SYSTEM
C
C   E1:  A(1,1) X(1) + A(1,2) X(2)                           = A(1,N+1)
C   E2:  A(2,1) X(1) + A(2,2) X(2) + A(2,3) X(3)             = A(2,N+1)
C   :
C   .
C   EN:                 A(N-1,N) X(N-1) + A(N,N) X(N)        = A(N,N+1)
C
C
C     INPUT:   THE DIMENSION N; THE ENTRIES OF A.
C
C     OUTPUT:  THE SOLUTION X(1),..,X(N).
C
C     INITIALIZATION
      DIMENSION A(10),B(10),C(10),BB(10),Z(10),X(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is Crout Method for tridiagonal'
      WRITE(6,*) 'linear systems.'
      WRITE(6,*) 'The array will be input from a text file in the'
      WRITE(6,*) 'order:  all diagonal entries, all lower sub-'
      WRITE(6,*) 'diagonal entries, all upper sub-diagonal entries,'
      WRITE(6,*) 'inhomogeneous term.'
      WRITE(6,*) 'Place each group of entries on a line,'
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
         WRITE(6,*) 'Input the number of equations - an integer '
         WRITE(6,*)
         READ(5,*) N
         IF (N .GT. 0) THEN
            NN=N-1
            READ(INP,*) (A(I),I=1,N)
            READ(INP,*) (B(I),I=2,N)
            READ(INP,*) (C(I),I=1,NN)
            READ(INP,*) (BB(I),I=1,N)
            OK = .TRUE.
            CLOSE(UNIT=INP)
         ELSE
            WRITE(6,*) 'The number must be a positive integer'
         ENDIF
         GOTO 19
      ELSE
         WRITE(6,*) 'The program will end so the input file can '
         WRITE(6,*) 'be created. '
         OK = .FALSE.
      ENDIF
111   IF(.NOT. OK) GOTO 400
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
      WRITE(OUP,*) 'CROUT METHOD FOR TRIDIAGONAL SYSTEMS'
      WRITE(OUP,4)
      WRITE(OUP,5) A(1),C(1),BB(1)
      WRITE(OUP,6) ( B(I),I-1,A(I),I,C(I),I+1,BB(I) ,I=2,NN)
      WRITE(OUP,7) B(N),A(N),BB(N)
C     STEP 1
C     THE ENTRIES OF U OVERWRITE C AND THE ENTRIES OF L OVERWRITE A
      C(1)=C(1)/A(1)
C     STEP 2
      DO 10 I=2,NN
C          ITH ROW OF L
           A(I)=A(I)-B(I)*C(I-1)
C          (I+1)ST COLUMN OF U
10         C(I)=C(I)/A(I)
C     STEP 3
C     NTH ROW OF L
      A(N)=A(N)-B(N)*C(N-1)
C     STEP 4
C     STEPS 4,5 SOLVE LZ = B
      Z(1)=BB(1)/A(1)
C     STEP 5
      DO 20 I=2,N
20         Z(I)=(BB(I)-B(I)*Z(I-1))/A(I)
C     STEP 6
C     STEPS 6,7 SOLVE UX = Z
      X(N)=Z(N)
C     STEP 7
      DO 30 II=1,NN
           I=NN-II+1
30         X(I)=Z(I)-C(I)*X(I+1)
C     STEP 8
      WRITE(OUP,8)
      WRITE(OUP,9) (X(I),I=1,N)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
4     FORMAT(1X,'ORIGINAL SYSTEM')
5     FORMAT(1X,E15.8,' X( 1) + ',E15.8,' X( 2) ',25X,'= ',E15.8)
6     FORMAT(1X,E15.8,' X(',I2,') + ',E15.8,' X(',I2,') + ',E15.8,' X('
     *,I2,') = ',E15.8)
7     FORMAT(1X,23X,E15.8,' X(N-1) + ',E15.8,' X( N) = ',E15.8)
8     FORMAT(1X,'THE SOLUTION VECTOR IS')
9     FORMAT(1X,4(5X,E15.8))
      END
