C***********************************************************************
C                                                                      *
C                FAST FOURIER TRANSFORM                                *
C                                                                      *
C***********************************************************************
C
C
C
C     TO COMPUTE THE COEFFICEINTS IN THE DISCRETE APPROXIMATION
C     F(X) FOR THE DATA (X(J),Y(J)), 0<=J<=2M-1 WHERE M=2**P AND
C     X(J) = -PI + J*PI/M  FOR 0 <= J <= 2M-1.
C
C     INPUT:   M; Y(0),Y(1),...,Y(2M-1).
C
C     OUTPUT:  COMPLEX NUMBERS C(0),...,C(2M-1); REAL NUMBERS
C              A(0), ..., A(M); B(1), ..., B(M-1)
C     NOTE:  THE MULTIPLICATION BY EXP(-K*PI*I) IS DONE WITHIN THE
C            THE PROGRAM.
C
      COMPLEX C(16), W(16), WW, T1, T3
      REAL Y(16)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      F(ZZ) = EXP(-1.0-ZZ/PI)
      PI = 3.1415926
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Fast Fourier Transform.'
      WRITE(6,*) ' '
      WRITE(6,*) 'The user must make provisions if the'
      WRITE(6,*) 'interval is not [-pi,pi].'
      WRITE(6,*) 'The example illustrates the required'
      WRITE(6,*) 'provisions under input method 3.'
      OK = .FALSE.
11    IF ( .NOT. OK) THEN
         WRITE(6,*) 'Choice of input method: '
         WRITE(6,*) '1. Input entry by entry from keyboard '
         WRITE(6,*) '2. Input data from a text file '
         WRITE(6,*) '3. Generate data using a function F'
         WRITE(6,*) 'Choose 1, 2, or 3 please '
         WRITE(6,*) ' '
         READ(5,*)  FLAG
         IF( ( FLAG .GE. 1 ) .AND. ( FLAG .LE. 3 )) OK = .TRUE.
         GOTO 11
      ENDIF
      IF (FLAG .EQ. 1) THEN
         OK = .FALSE.
21       IF (.NOT. OK ) THEN
            WRITE(6,*) 'Input number m '
            WRITE(6,*) ' '
            READ(5,*) M
            IF (M .GT. 0 ) THEN
               OK = .TRUE.
               N=2*M
               DO 31 I = 1, N
                  J=I-1
                  WRITE(6,*) 'Input Y(',J,') '
                  WRITE(6,*) ' '
                  READ(5,*) Y(I)
31             CONTINUE
            ELSE
               WRITE(6,*) 'Number must be a positive integer '
            ENDIF
            GOTO 21
         ENDIF
      ENDIF
      IF (FLAG .EQ. 2) THEN
         WRITE(6,*) 'Has a text file been created with the data'
         WRITE(6,*) 'y(0), ..., y(2m-1) separated by blanks? '
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
41          IF (.NOT. OK) THEN
               WRITE(6,*) 'Input number m '
               WRITE(6,*) ' '
               READ(5,*) M
               IF ( M .GT. 0) THEN
                  OK = .TRUE.
                  N=2*M
                  READ(4,*) (Y(I),I=1,N)
                  CLOSE(UNIT=4)
               ELSE
                  WRITE(6,*) 'Number must be a positive integer '
               ENDIF
               GOTO 41
            ENDIF
         ELSE
            WRITE(6,*) 'The program will end so the input file can '
            WRITE(6,*) 'be created. '
            OK = .FALSE.
         ENDIF
      ENDIF
      IF (FLAG .EQ. 3) THEN
         WRITE(6,*) 'Has the function F been created in the program '
         WRITE(6,*) 'Enter Y or N - letter within quotes'
         WRITE(6,*) ' '
         READ(5,*)  AA
         IF (( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
            OK = .FALSE.
61          IF (.NOT. OK) THEN
               WRITE(6,*) 'Input number m '
               WRITE(6,*) ' '
               READ(5,*) M
               IF ( M .GT. 0 ) THEN
                  OK = .TRUE.
                  N=2*M
                  DO 71 I = 1, N
                     XX=(I-1.0)/8.0
                     Z = PI*(XX-1)
                     Y(I) = F(Z)
71                CONTINUE
               ELSE
                  WRITE(6,*) 'Number must be a positive integer '
               ENDIF
               GOTO 61
            ENDIF
         ELSE
            WRITE(6,*) 'The program will end so that the function F '
            WRITE(6,*) 'can be created '
            OK = .FALSE.
         ENDIF
      ENDIF
      IF (.NOT. OK) GOTO 500
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
      WRITE(OUP,*) 'FAST FOURIER TRANSFORM'
C     STEP 1
C     USE N2 FOR M, NG FOR P, NU1 FOR Q, WW FOR ZETA
      N2 = N/2
      NG = ALOG(FLOAT(N))/ALOG(FLOAT(2)) + .5
      NU1 = NG - 1
      WW = CEXP ( CMPLX( 0.0 , 2*PI/N))
C     STEP 2
C     SUBSCRIPTS ARE SHIFTED TO AVOID ZERO SUBSCRIPTS
      DO 50 I=1,N
50         C(I) = CMPLX( Y(I), 0.0)
C     STEP 3
      DO 40 I=1,N2
           W(I)=WW**I
40         W(N2+I)=-W(I)
C     STEP 4
      K=0
C     STEP 5
      DO 20 L=1,NG
C          STEP 6
100        IF (K .GE. N-1) GOTO 200
C               STEP 7
                DO 30 I=1,N2
C                    STEP 8
                     M1=K/2**NU1
C                    THE SUBPROGRAM IBR DOES THE BIT REVERSAL
                     NP=IBR(M1,NG)
C                    T1 IS THE SAME AS ETA
                     T1=C(K+N2+1)
C                    STEP 9
                     IF(NP.NE.0) T1=T1*W(NP)
                     C(K+N2+1)=C(K+1)-T1
                     C(K+1)=C(K+1)+T1
C                    STEP 10
30                   K=K+1
C               STEP 11
                K=K+N2
           GOTO 100
C          STEP 12
200        K=0
           N2=N2/2
20         NU1=NU1-1
C     STEP 13
300   IF (K .GE. N-1) GOTO 400
C          STEP 14
           I=IBR(K,NG)
C          STEP 15
           IF(I.GT.K) THEN
                T3=C(K+1)
                C(K+1)=C(I+1)
                C(I+1)=T3
           END IF
C          STEP 16
           K=K+1
      GOTO 300
C     STEP 17 AND 18
400   DO 60 I=1,N
           C(I)=CEXP(CMPLX(0.0,-(I-1)*PI))*C(I)/(FLOAT(N)/2)
60         WRITE(OUP,2) I,C(I)
500   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
2     FORMAT(1X,I2,2(2X,E15.8))
      END
C
           FUNCTION IBR(J,NU)
           J1=J
           IBR=0
           DO 70 I=1,NU
                J2=J1/2
                IBR=IBR*2+(J1-2*J2)
70         J1=J2
           RETURN
           END
