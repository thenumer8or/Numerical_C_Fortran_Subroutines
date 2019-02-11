C***********************************************************************
C                                                                      *
C               POISSON FINITE-DIFFERENCE METHOD                       *
C                                                                      *
C***********************************************************************
C
C
C
C     TO APPROXIMATE THE SOLUTION TO THE POISSON EQUATION
C         DEL(U) = F(X,Y), A <= X <= B, C <= Y <= D,
C     SUBJECT TO BOUNDARY CONDITIONS:
C                     U(X,Y)=G(X,Y),
C     IF X = A OR X = B, FOR C <= Y <= D
C     IF Y = C OR Y = D, FOR A <= X <= B
C
C     INPUT:   ENDPOINTS A, B, C, D; INTEGERS M, N; TOLERENCE TOL;
C              MAXIMUM NUMBER OF ITERATIONS M.
C
C     OUTPUT:  APPROXIMATIONS W(I,J) TO U(X(I),Y(J)) FOR EACH
C              I=1,...,N-1 AND J=1,...,M-1 OR A MESSAGE THAT THE
C              MAXIMUM NUMBER OF ITERATIONS WAS EXCEEDED.
C
C     INITIALIZATION
C     DIMENSION W(N-1,M-1),X(N-1),Y(M-1)
      DIMENSION W(5,4),X(5),Y(4)
      CHARACTER NAME*14,NAME1*14,AA*1
      INTEGER INP,OUP,FLAG,LBOUND
      LOGICAL OK
C     USE  -F(X,Y) INSTEAD OF F(X,Y) AND WE USED EXACT SOLUTION FOR
C     BOUNDARY VALUES
      F(XZ,YZ)=-XZ*EXP(YZ)
      G(XZ,YZ)= XZ*EXP(YZ)
      OPEN(UNIT=5,FILE='CON',ACCESS='SEQUENTIAL')
      OPEN(UNIT=6,FILE='CON',ACCESS='SEQUENTIAL')
      WRITE(6,*) 'This is the Linear Finite Difference Method'
      WRITE(6,*) 'for Elliptic Equations.'
      WRITE(6,*) 'Have the functions F(x,y) and G(x,y) been created?'
      WRITE(6,*) 'Enter Y or N '
      WRITE(6,*) ' '
      READ(5,*)  AA
      IF(( AA .EQ. 'Y' ) .OR. ( AA .EQ. 'y' )) THEN
         OK = .FALSE.
19       IF (OK) GOTO 11
            WRITE(6,*) 'Input endpoints [A,B] on X-axis'
            WRITE(6,*) 'separated by blank.'
            WRITE(6,*) ' '
            READ(5,*) A, B
            WRITE(6,*) 'Input endpoints [C,D] on X-axis'
            WRITE(6,*) 'separated by blank.'
            WRITE(6,*) ' '
            READ(5,*) C, D
            IF ((A.GE.B).OR.(C.GE.D)) THEN
               WRITE(6,*) 'Left endpoint must be less'
               WRITE(6,*) 'than right endpoint'
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 19
11       OK = .FALSE.
12       IF (OK) GOTO 13
            WRITE(6,*) 'Input number of intervals n on the X-axis'
            WRITE(6,*) 'and m on the Y-axis separated by a blank.'
            WRITE(6,*) 'Note:  both n and m should be larger than 2.'
            WRITE(6,*) ' '
            READ(5,*) N, M
            IF ((N.LE.0).OR.(M.LE.0)) THEN
              WRITE(6,*) 'Must be positive integers '
            ELSE
              OK = .TRUE.
            ENDIF
         GOTO 12
13       OK = .FALSE.
14       IF (OK) GOTO 15
            WRITE(6,*) 'Input tolerance '
            WRITE(6,*) ' '
            READ(5,*) TOL
            IF (TOL.LE.0.0) THEN
               WRITE(6,*) 'Tolerance must be positive '
            ELSE
               OK = .TRUE.
            ENDIF
         GOTO 14
15       OK = .FALSE.
16       IF(OK) GOTO 17
            WRITE(6,*) 'Input the maximum number of iterations.'
            WRITE(6,*) ' '
            READ(5,*) LBOUND
            IF(LBOUND.LE.0) THEN
               WRITE(6,*) 'Must be a positive integer.'
            ELSE
               OK = .TRUE.
            ENDIF
            GO TO 16
17       CONTINUE
      ELSE
         WRITE(6,*) 'The program will end so that the functions'
         WRITE(6,*) 'F and G can be created.'
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
      WRITE(OUP,*) 'POISSON EQUATION FINITE-DIFFERENCE METHOD'
      WRITE(OUP,*) A,B,C,D
      WRITE(OUP,*) TOL,M,N,LBOUND
      MM=M-1
      MMM=M-2
      NN=N-1
      NNN=N-2
C     XK WILL BE USED FOR K
C     STEP 1
      H=(B-A)/N
      XK=(D-C)/M
C     A,B,C,D WILL BE USED FOR X(0),X(N),Y(0),Y(M)
C     STEPS 2, 3 CONSTRUCT MESH POINTS
C     STEP 2
      DO 10 I=1,NN
10         X(I)=A+I*H
C     STEP 3
      DO 20 J=1,MM
20         Y(J)=C+J*XK
C     STEP 4
      DO 30 I=1,NN
           DO 30 J=1,MM
30              W(I,J)=0
C     STEP 5
C     USE V FOR LAMBDA, VV FOR MU
      V=H*H/(XK*XK)
      VV=2*(1+V)
      L=1
C     Z IS A NEW VALUE OF W(I,J) TO BE USED IN COMPUTING THE ERROR E
C     USE E INSTEAD OF NORM
C     STEP 6
100   IF (L.GT.LBOUND) GOTO 200
C     STEPS 7 THROUGH 20 PERFORM GAUSS-SEIDEL ITERATIONS
C          STEP 7
           Z=(H*H*F(X(1),Y(MM))+G(A,Y(MM))+V*G(X(1),D)+W(1,MM-1)*V+W(2,
     *     MM))/VV
           E=ABS(W(1,MM)-Z)
           W(1,MM)=Z
C          STEP 8
           DO 40 I=2,NNN
                Z=(H*H*F(X(I),Y(MM))+V*G(X(I),D)+W(I-1,MM)+W(I+1,MM)+V*
     *          W(I,MM-1))/VV
                IF(ABS(W(I,MM)-Z).GT.E) E=ABS(W(I,MM)-Z)
40         W(I,MM)=Z
C          STEP 9
           Z=(H*H*F(X(NN),Y(MM))+G(B,Y(MM))+V*G(X(NN),D)+W(NN-1,MM)+V*
     *     W(NN,MM-1))/VV
           IF(ABS(W(NN,MM)-Z).GT.E) E=ABS(W(NN,MM)-Z)
           W(NN,MM)=Z
C          STEP 10
           DO 50 LL=2,MMM
                J=MMM-LL+2
C               STEP 11
                Z=(H*H*F(X(1),Y(J))+G(A,Y(J))+V*W(1,J+1)+V*W(1,J-1)+W(2
     *          ,J))/VV
                IF(ABS(W(1,J)-Z).GT.E) E=ABS(W(1,J)-Z)
                W(1,J)=Z
C               STEP 12
                DO 60 I=2,NNN
                     Z=(H*H*F(X(I),Y(J))+W(I-1,J)+V*W(I,J+1)+V*W(I,J-1)
     *               +W(I+1,J))/VV
                     IF(ABS(W(I,J)-Z).GT.E) E=ABS(W(I,J)-Z)
60              W(I,J)=Z
C               STEP 13
                Z=(H*H*F(X(NN),Y(J))+G(B,Y(J))+W(NN-1,J)+V*W(NN,J+1)+V*
     *          W(NN,J-1))/VV
                IF(ABS(W(NN,J)-Z).GT.E) E=ABS(W(NN,J)-Z)
50         W(NN,J)=Z
C          STEP 14
           Z=(H*H*F(X(1),Y(1))+V*G(X(1),C)+G(A,Y(1))+V*W(1,2)+W(2,1))/VV
           IF(ABS(W(1,1)-Z).GT.E) E=ABS(W(1,1)-Z)
           W(1,1)=Z
C          STEP 15
           DO 70 I=2,NNN
                Z=(H*H*F(X(I),Y(1))+V*G(X(I),C)+W(I+1,1)+W(I-1,1)+V*W(I
     *          ,2))/VV
                IF(ABS(W(I,1)-Z).GT.E) E=ABS(W(I,1)-Z)
70         W(I,1)=Z
C          STEP 16
           Z=(H*H*F(X(NN),Y(1))+V*G(X(NN),C)+G(B,Y(1))+W(NN-1,1)+V*W(NN
     *     ,2))/VV
           IF(ABS(W(NN,1)-Z).GT.E) E=ABS(W(NN,1)-Z)
           W(NN,1)=Z
C          STEP 17
           IF(E.LE.TOL) THEN
                WRITE(OUP,2) L
C               STEP 18
                DO 80 I=1,NN
                DO 80 J=1,MM
80                   WRITE(OUP,3) I,J,X(I),Y(J),W(I,J)
C               STEP 19
                GOTO 400
           END IF
C          STEP 20
           L=L+1
      GOTO 100
C     STEP 21
200   WRITE(OUP,1) LBOUND
400   CLOSE(UNIT=5)
      CLOSE(UNIT=OUP)
      IF(OUP.NE.6) CLOSE(UNIT=6)
      STOP
1     FORMAT(1X,'METHOD FAILS AFTER ITERATION NO. ',I4)
2     FORMAT(1X,'NO. OF ITERATIONS IS',I4,3X,'ORDER OF OUTPUT IS I,J,X(I
     *),Y(J),W(I,J)')
3     FORMAT(1X,2(I3,2X),3(E15.8,3X))
      END
