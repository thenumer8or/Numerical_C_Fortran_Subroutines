C***********************************************************************
C                                                                      *
C              NONLINEAR FINITE-DIFFERENCE METHOD                      *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION TO THE NONLINEAR BOUNDARY-VALUE
C     PROBLEM
C
C         Y'' = F(X,Y,Y'), A <= X <= B, Y(A) = ALPHA, Y(B) = BETA:
C
C     INPUT:   ENDPOINTS A,B; BOUNDARY CONDITIONS ALPHA, BETA;
C              INTEGER N; TOLERANCE TOL; MAXIMUM NUMBER OF ITERATIONS
C              M.
C
C     OUTPUT:  APPROXIMATIONS W(I) TO Y(X(I)) FOR EACH I=0,1,...N+1
C              OR A MESSAGE THAT THE MAXIMUM NUMBER OF ITERATIONS WAS
C              EXCEEDED.
C
C     INITIALIZATION
      DIMENSION W(20),A(20),C(20),D(20),XL(20),XU(20),Z(20),V(20),B(20)
      CHARACTER NAME*14,NAME1*14,AAA*1
      INTEGER INP,OUP,FLAG
      LOGICAL OK
      F(XZ,YZ,ZZ)=(32+2*XZ**3-YZ*ZZ)/8
C     FY(XZ,YZ,ZZ) REPRESENTS PARTIAL OF F WITH RESPECT TO Y
      FY(XZ,YZ,ZZ)=-ZZ/8
C     FYP(XZ,YZ,ZZ) REPRESENTS PARTIAL OF F WITH RESPECT TO Y'
      FYP(XZ,YZ,ZZ)=-YZ/8
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Nonlinear Finite Difference Method.'
      WRITE(6,*) 'Have the functions F, FY (partial of F'
      WRITE(6,*) 'with respect to y) and FYP (partial of'
      WRITE(6,*) 'F with respect to y'') been created? '
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AAA
      IF(( AAA .EQ. 'Y' ) .OR. ( AAA .EQ. 'y' )) THEN
         OK = .FALSE.
110      IF (OK) GOTO 11
            WRITE(6,*) 'Input left and right endpoints separated by'
            WRITE(6,*) 'blank - include decimal point'
            WRITE(6,*) ' '
            READ(5,*) AA, BB
            IF (AA.GE.BB) THEN
               WRITE(6,*) 'Left endpoint must be less'
               WRITE(6,*) 'than right endpoint'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 110
11       OK = .FALSE.
         WRITE(6,*) 'Input Y(',AA,')'
         WRITE(6,*) ' '
         READ(5,*) ALPHA
         WRITE(6,*) 'Input Y(',BB,')'
         WRITE(6,*) ' '
         READ(5,*) BETA
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input a positive integer for the number'
            WRITE(6,*) 'of subinvervals.  note: h = (b-a)/(n+1)'
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
      WRITE(OUP,*) 'NONLINEAR FINITE DIFFERENCE METHOD'
C     STEP 1
      N1=N-1
      H=(BB-AA)/(N+1)
C     STEP 2
      DO 10 I=1,N
10    W(I)=ALPHA + I*H*(BETA-ALPHA)/(BB-AA)
C     STEP 3
      K=1
C     STEP 4
100   IF (K.GT.M) GOTO 200
C          STEP 5
           X=AA+H
           T=(W(2)-ALPHA)/(2*H)
           A(1)=2+H*H*FY(X,W(1),T)
           B(1)=-1+H*FYP(X,W(1),T)/2
           D(1)=-(2*W(1)-W(2)-ALPHA+H*H*F(X,W(1),T))
C          STEP 6
           DO 20 I=2,N1
                X=AA+I*H
                T=(W(I+1)-W(I-1))/(2*H)
                A(I)=2+H*H*FY(X,W(I),T)
                B(I)=-1+H*FYP(X,W(I),T)/2
                C(I)=-1-H*FYP(X,W(I),T)/2
20         D(I)=-(2*W(I)-W(I+1)-W(I-1)+H*H*F(X,W(I),T))
C          STEP 7
           X=BB-H
           T=(BETA-W(N-1))/(2*H)
           A(N)=2+H*H*FY(X,W(N),T)
           C(N)=-1-H*FYP(X,W(N),T)/2
           D(N)=-(2*W(N)-W(N-1)-BETA+H*H*F(X,W(N),T))
C          STEP 8
C          STEPS 8 THROUGH 14 SOLVE A TRIADIAGONAL LINEAR SYSTEM
C          USING ALGORITHM 6.7
C          USE XL FOR L AND XU FOR U
           XL(1)=A(1)
           XU(1)=B(1)/A(1)
C          STEP 9
           DO 30 I=2,N1
                XL(I)=A(I)-C(I)*XU(I-1)
30         XU(I)=B(I)/XL(I)
C          STEP 10
           XL(N)=A(N)-C(N)*XU(N-1)
C          STEP 11
           Z(1)=D(1)/XL(1)
C          STEP 12
           DO 40 I=2,N
40         Z(I)=(D(I)-C(I)*Z(I-1))/XL(I)
C          STEP 13
           V(N)=Z(N)
C          VMAX IS USED TO MEASURE THE INFINITY NORM OF V
           VMAX=ABS(V(N))
           W(N)=W(N)+V(N)
C          STEP 14
           DO 50 J=1,N1
                I=N-J
                V(I)=Z(I)-XU(I)*V(I+1)
                W(I)=W(I)+V(I)
                IF(ABS(V(I)).GT.VMAX) VMAX=ABS(V(I))
50         CONTINUE
C          STEP 15
C          TEST FOR ACCURACY
           IF(VMAX.LE.TOL) THEN
C               STEP 16
                WRITE(OUP,2) K
                WRITE(OUP,3) AA,ALPHA
                DO 60 I=1,N
                     X=AA+I*H
60              WRITE(OUP,3) X,W(I)
                WRITE(OUP,3) BB,BETA
C               STEP 17
C               PROCEDURE COMPLETED SUCCESSFULLY
                GOTO 400
           END IF
C          STEP 18
           K=K+1
      GOTO 100
C     STEP 19
200   WRITE(OUP,1) M
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'FAILURE AFTER ITERATION NO.',I4)
2     FORMAT(1X,'ORDER OF OUTPUT - X(I),W(I)',1X,I3,1X,'ITERATIONS')
3     FORMAT(1X,2(E15.8,3X))
      END
