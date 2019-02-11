C***********************************************************************
C                                                                      *
C                   NONLINEAR SHOOTING METHOD                          *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION OF THE NONLINEAR BOUNDARY-VALUE
C     PROBLEM
C
C            Y'' = F(X,Y,Y'), A<=X<=B, Y(A)=ALPHA, Y(B)=BETA:
C
C
C     INPUT:   ENDPOINTS A,B; BOUNDARY CONDITIONS ALPHA, BETA; NUMBER OF
C              SUBINTERVALS N; TOLERANCE TOL; MAXIMUM NUMBER OF
C              ITERATIONS M.
C
C     OUTPUT:  APPROXIMATIONS W(1,I) TO Y(X(I)); W(2,I) TO Y'(X(I))
C              FOR EACH I=0,1,...,N OR A MESSAGE THAT THE MAXIMUM
C              NUMBER OF ITERATIONS WAS EXCEEDED.
C
C     INITIALIZATION
      DIMENSION W1(21),W2(21)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      F(X,Y,Z)=(32+2*X**3-Y*Z)/8
C     FY(X,Y,Z) REPRESENTS PARTIAL OF F WITH RESPECT TO Y
      FY(X,Y,Z)=-Z/8
C     FYP(X,Y,Z) REPRESENTS PARTIAL OF F WITH RESPECT TO Y'
      FYP(X,Y,Z)=-Y/8
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Nonlinear Shooting Method.'
      WRITE(6,*) 'Have the functions F, FY (partial of F'
      WRITE(6,*) 'with respect to y) and FYP (partial of'
      WRITE(6,*) 'F with respect to y'') been created? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
110      IF (OK) GOTO 11
            WRITE(6,*) 'Input left and right endpoints separated by'
            WRITE(6,*) 'blank - include decimal point'
            WRITE(6,*) ' '
            READ(5,*) A, B
            IF (A.GE.B) THEN
               WRITE(6,*) 'Left endpoint must be less'
               WRITE(6,*) 'than right endpoint'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 110
11       OK = .FALSE.
         WRITE(6,*) 'Input Y(',A,')'
         WRITE(6,*) ' '
         READ(5,*) ALPHA
         WRITE(6,*) 'Input Y(',B,')'
         WRITE(6,*) ' '
         READ(5,*) BETA
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input a positive integer for the number'
            WRITE(6,*) 'of subinvervals '
            WRITE(6,*) ' '
            READ(5,*) N
            IF ( N .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       OK=.FALSE.
14       IF (OK) GOTO 15
            WRITE(6,*) 'Input tolerance '
            WRITE(6,*) ' '
            READ(5,*) TOL
            IF (TOL.LE.0.0) THEN
               WRITE(6,*) 'Tolerance must be positive '
               WRITE(6,*) ' '
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 14
15       OK = .FALSE.
16       IF (OK) GOTO 17
            WRITE(6,*) 'Input maximum number of iterations '
            WRITE(6,*) '- no decimal point '
            WRITE(6,*) ' '
            READ(5,*) M
            IF ( M .LE. 0 ) THEN
              WRITE(6,*) 'Must be positive integer '
              WRITE(6,*) ' '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 16
17       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'F, FY and FYP can be created '
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
      WRITE(OUP,*) 'NONLINEAR SHOOTING METHOD'
C     STEP 1
      H=(B-A)/N
      K=1
      TK=(BETA-ALPHA)/(B-A)
C     STEP 2
100   IF (K.GT.M) GOTO 200
C          STEP 3
C          W1(1)=W(1,0), W2(1)=W(2,0), U1=U(1), U2=U(2)
           W1(1)=ALPHA
           W2(1)= TK
           U1=0
           U2=1
C          STEP 4
C          RUNGE-KUTTA METHOD FOR SYSTEMS IS USED IN STEPS 5, 6
           DO 10 I=1,N
C               STEP 5
                X=A+(I-1)*H
                T=X+H/2
C               STEP 6
                XK11=H*W2(I)
                XK12=H*F(X,W1(I),W2(I))
                XK21=H*(W2(I)+XK12/2)
                XK22=H*F(T,W1(I)+XK11/2,W2(I)+XK12/2)
                XK31=H*(W2(I)+XK22/2)
                XK32=H*F(T,W1(I)+XK21/2,W2(I)+XK22/2)
                XK41=H*(W2(I)+XK32)
                XK42=H*F(X+H,W1(I)+XK31,W2(I)+XK32)
                W1(I+1)=W1(I)+(XK11+2*XK21+2*XK31+XK41)/6
                W2(I+1)=W2(I)+(XK12+2*XK22+2*XK32+XK42)/6
                YK11=H*U2
                YK12=H*(FY(X,W1(I),W2(I))*U1+FYP(X,W1(I),W2(I))*U2)
                YK21=H*(U2+YK12/2)
                YK22=H*(FY(T,W1(I),W2(I))*(U1+YK11/2)
     *          +FYP(T,W1(I),W2(I))*(U2+YK21/2))
                YK31=H*(U2+YK22/2)
                YK32=H*(FY(T,W1(I),W2(I))*(U1+YK21/2)
     *          +FYP(T,W1(I),W2(I))*(U2+YK22/2))
                YK41=H*(U2+YK32)
                YK42=H*(FY(X+H,W1(I),W2(I))*(U1+YK31)+FYP(X+H,W1(I),
     *          W2(I))*(U2+YK32))
                U1=U1+(YK11+2*YK21+2*YK31+YK41)/6
                U2=U2+(YK21+2*YK22+2*YK32+YK42)/6
10         CONTINUE
C          TEST FOR ACCURACY
C          STEP 7
           IF(ABS(W1(N+1)-BETA).LE.TOL) THEN
C               PROCESS IS COMPLETE
C               W1(J), W2(J) ARE THE APPROXIMATIONS FOR Y AND Y'
                WRITE(OUP,2) K
C               STEP 8
                DO 20 I=1,N
                     J=I+1
                     X=A+I*H
20              WRITE(OUP,1) X,W1(J),W2(J)
C               STEP 9
                GOTO 400
           END IF
C          STEP 10
C          NEWTON'S METHOD APPLIED TO IMPROVE TK
           TK=TK-(W1(N+1)-BETA)/U1
           K=K+1
      GOTO 100
C     STEP 11
C     METHOD FAILED
200   WRITE(OUP,3) M
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,3(E15.8,3X))
2     FORMAT(1X,'ORDER OF OUTPUT - X(I),W1(I),W2(I)',1X,I3,1X,'ITERATION
     *S')
3     FORMAT(1X,'METHOD FAILED AFTER ITERATION NO.' ,I4)
      END
