C***********************************************************************
C                                                                      *
C                   CHOLESKI'S  METHOD                                 *
C                                                                      *
C***********************************************************************
C
C
C
C     TO FACTOR THE POSITIVE DEFINITE N BY N MATRIX A INTO LL**T,
C     WHERE L IS LOWER TRIANGULAR.
C
C     INPUT:   THE DIMENSION N; ENTRIES A(I,J), 1<=I, J<=N OF A.
C
C     OUTPUT:  THE ENTRIES L(I,J), 1<=J<=I, 1<=I<=N OF L.
C
C     THE ENTRIES OF U=L**T ARE U(I,J)=L(J,I), I<=J<=N, 1<=I<=N
C
C     INITIALIZATION
      DIMENSION A(10,10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is Choleski Factorization Method'
      WRITE(6,*) 'The array will be input from a text file in the'
      WRITE(6,*) ' order: A(1,1), A(1,2), ..., A(1,N), A(2,1),'
      WRITE(6,*) ' A(2,2), ..., A(2,N)..., A(N,1), A(N,2),'
      WRITE(6,*) ' ..., A(N,N) '
      WRITE(6,*) 'Place as many entries as desired on each line,'
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
         WRITE(6,*) 'Input the number of rows - an integer '
         WRITE(6,*)
         READ(5,*) N
         IF (N .GT. 0) THEN
            READ(INP,*) ((A(I,J), J=1,N),I=1,N)
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
      WRITE(OUP,*) 'CHOLESKI FACTORIZATION'
      WRITE(OUP,4)
      WRITE(OUP,5) ((A(I,J),J=1,N),I=1,N)
C     STEP 1
      A(1,1)=SQRT(A(1,1))
C     STEP 2
      DO 10 J=2,N
10    A(J,1)=A(J,1)/A(1,1)
C     STEP 3
      NN=N-1
      DO 20 I=2,NN
C          STEP 4
           KK=I-1
           S=0
C          DO-LOOP COMPUTES THE SUMMATION
           DO 30 K=1,KK
30              S=S-A(I,K)**2
           A(I,I)=SQRT(A(I,I)+S)
C          STEP 5
           JJ=I+1
           DO 40 J=JJ,N
                S=0
                KK=I-1
C               DO-LOOP COMPUTES THE SUMMATION
                DO 50 K=1,KK
50                   S=S-A(J,K)*A(I,K)
40              A(J,I)=(A(J,I)+S)/A(I,I)
20         CONTINUE
C     STEP 6
      S=0
      DO 60 K=1,NN
60         S=S-A(N,K)**2
      A(N,N)=SQRT(A(N,N)+S)
C     STEP 7
      WRITE(OUP,6)
      DO 70 I=1,N
70    WRITE(OUP,5) (A(I,J),J=1,I)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
4     FORMAT(1X,'ORIGINAL MATRIX A')
5     FORMAT(1X,3(5X,E15.8))
6     FORMAT(1X,'MATRIX L')
      END
