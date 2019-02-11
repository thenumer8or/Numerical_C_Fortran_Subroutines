C***********************************************************************
C                                                                      *
C                    LDL-T  FACTORIZATION                              *
C                                                                      *
C***********************************************************************
C
C
C
C     TO FACTOR THE POSITIVE DEFINITE N BY N MATRIX A INTO LDL**T,
C     WHERE L IS LOWER TRIANGULAR WITH ONES ALONG THE DIAGONAL AND D
C     IS A DIAGONAL MATRIX WITH POSITIVE ENTRIES ON THE DIAGONAL
C
C     INPUT:   THE DIMENSION N; ENTRIES A(I,J), 1<=I, J<=N OF A.
C
C     OUTPUT:  THE ENTRIES L(I,J), 1<=J<I, 1<=I<=N OF L AND D(I),
C              1<=I<=N OF D.
C
C
C
C     INITIALIZATION
      DIMENSION A(10,10),D(10),V(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the LDL^t method for Positive Definite'
      WRITE(6,*) 'Matrices.'
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
      WRITE(OUP,*) 'LDL^t METHOD FOR POSITIVE DEFINITE MATRICES'
      WRITE(OUP,4)
      WRITE(OUP,5) ((A(I,J),J=1,N),I=1,N)
C     STEP 1
      DO 50 I = 1, N
C         STEP 2
          II = I-1
          IF (I .GT. 1 ) THEN
                DO 10 J=1,II
10              V(J) = A(I,J)*D(J)
          ENDIF
C         STEP 3
          D(I) = A(I,I)
          IF (I .GT. 1 ) THEN
                DO 20 J=1,II
20              D(I) = D(I)-A(I,J)*V(J)
          ENDIF
C         STEP 4
          JJ = I+1
          IF (I .LT. N) THEN
                DO 30 J=JJ,N
                     IF (I .GT. 1) THEN
                          DO 40 K=1,II
40                        A(J,I)=A(J,I)-A(J,K)*V(K)
                     ENDIF
30              A(J,I)=A(J,I)/D(I)
          ENDIF
50    CONTINUE
C     STEP 5
      WRITE(OUP,6)
      DO 70 I=2,N
          II=I-1
70    WRITE(OUP,5) (A(I,J),J=1,II)
      WRITE(OUP,7)
      WRITE(OUP,5) (D(I),I=1,N)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
4     FORMAT(1X,'ORIGINAL MATRIX A')
5     FORMAT(1X,3(5X,E15.8))
6     FORMAT(1X,'MATRIX L')
7     FORMAT(1X,'DIAGONAL D')
      END
