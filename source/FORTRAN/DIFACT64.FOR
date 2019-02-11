C***********************************************************************
C                                                                      *
C              DIRECT FACTORIZATION                                    *
C                                                                      *
C***********************************************************************
C
C
C
C     TO FACTOR THE N BY N MATRIX A=(A(I,J)) INTO THE PRODUCT OF THE
C     LOWER TRIANGULAR MATRIX L= L(I,J)  AND U= U(I,J) , THAT IS A=LU,
C     WHERE THE MAIN DIAGONAL OF EITHER L OR U IS GIVEN:
C
C     INPUT:   DIMENSION N; THE ENTRIES A(I,J), 1<=I, J<=N, OF A; THE
C              DIAGONAL L(1,1),...L(N,N) OF L OR THE DIAGONAL
C              U(1,1),...U(N,N) OF U.
C
C     OUTPUT:  THE ENTRIES L(I,J), 1<=J<=I, 1<=I<=N OF L AND THE
C              ENTRIES U(I,J), I<=J<=N, 1<=I<=N OF U.
C
C     INITIALIZATION
      DIMENSION A(10,10),XL(10)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the LU factorization method'
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
            WRITE(6,*) 'Choice of diagonals:'
            WRITE(6,*) '1. Diagonal of L consists of ones'
            WRITE(6,*) '2. Diagonal of U consists of ones'
            WRITE(6,*) 'Please enter 1 or 2'
            WRITE(6,*) ' '
            READ(5,*) FLAG
            ISW=1
            IF(FLAG.EQ.1) ISW = 0
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
      WRITE(OUP,*) 'GENERAL LU FACTORIZATION '
      WRITE(OUP,3)
      DO 80 I=1,N
80    XL(I) = 1.0
      WRITE(OUP,5)((A(I,J),J=1,N),I=1,N)
C     STEP 1
      IF ( ABS(A(1,1)) .LT. 1.0E-20 ) THEN
           WRITE(OUP,8)
           GOTO 400
      END IF
C     THE ENTRIES OF L BELOW THE MAIN DIAGONAL WILL BE PLACED IN THE
C     CORRESPONDING ENTRIES OF A; THE ENTRIES OF U ABOVE THE MAIN
C     DIAGONAL WILL BE PLACED IN THE CORRESPONDING ENTRIES OF A;
C     THE MAIN DIAGONAL WHICH WAS NOT INPUT WILL BECOME THE MAIN
C     DIAGONAL OF A; THE INPUT MAIN DIAGONAL OF L OR U IS, OF COURSE,
C     PLACED IN XL.
      A(1,1)=A(1,1)/XL(1)
C     STEP 2
      DO 10 J=2,N
      IF(ISW.EQ.0) THEN
C          FIRST ROW OF U
           A(1,J)=A(1,J)/XL(1)
C          FIRST COLUMN OF L
           A(J,1)=A(J,1)/A(1,1)
      ELSE
C          FIRST ROW OF U
           A(1,J)=A(1,J)/A(1,1)
C          FIRST COLUMN OF L
           A(J,1)=A(J,1)/XL(1)
      END IF
10    CONTINUE
C     STEP 3
      M=N-1
      DO 20 I=2,M
C          STEP 4
           KK=I-1
           S=0
           DO 30 K=1,KK
30         S=S-A(I,K)*A(K,I)
           A(I,I)=(A(I,I)+S)/XL(I)
           IF(ABS(A(I,I)).LT.1.0E-20) THEN
                WRITE(OUP,8)
                GOTO 400
           END IF
C          STEP 5
           JJ=I+1
           DO 40 J=JJ,N
                SS=0
                S=0
                DO 50 K=1,KK
                     SS=SS-A(I,K)*A(K,J)
50                   S=S-A(J,K)*A(K,I)
                IF(ISW.EQ.0) THEN
C                    ITH ROW OF U
                     A(I,J)=(A(I,J)+SS)/XL(I)
C                    ITH COLUMN OF L
                     A(J,I)=(A(J,I)+S)/A(I,I)
                ELSE
C                    ITH ROW OF U
                     A(I,J)=(A(I,J)+SS)/A(I,I)
C                    ITH COLUMN OF L
                     A(J,I)=(A(J,I)+S)/XL(I)
                END IF
40         CONTINUE
20    CONTINUE
C     STEP 6
      S=0
      DO 60 K=1,M
60    S=S-A(N,K)*A(K,N)
      A(N,N)=(A(N,N)+S)/XL(N)
C     IF A(N,N) = 0 THEN A = LU BUT THE MATRIX IS SINGULAR
C     PROCESS IS COMPLETE, ALL ENTRIES OF A HAVE BEEN DETERMINED
C     STEP 7
      WRITE(OUP,6)
      WRITE(OUP,11) ISW
      WRITE(OUP,5) (XL(I),I=1,N)
      WRITE(OUP,7)
      WRITE(OUP,5) ((A(I,J),J=1,N),I=1,N)
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
3     FORMAT(1X,'ORIGINAL MATRIX A')
5     FORMAT(4(5X,E15.8))
6     FORMAT(1X,'DIAGONAL OF L OR U')
7     FORMAT(1X,'ENTRIES OF L BELOW/ON DIAGONAL AND ENTRIES OF U ABOVE/
     *ON DIAGONAL')
8     FORMAT(1X,' FACTORIZATION IMPOSSIBLE')
11    FORMAT(1X,'IF 0 THEN L AND IF 1 THEN U, IT IS ',I2)
      END
